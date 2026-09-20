/*
 * gp_vehicle_model.c
 */

#include "gp_vehicle_model.h"

void gp_estimate_fz(float vx, float ax, float ay, float fz_out[4]) {
    const float total_weight = GP_MASS * GP_GRAVITY;
    const float inv_wb = 1.0f / (GP_WB + 1e-3f);

    // Reparto estático real por posición del CG (45% Delante / 55% Detrás)
    float fz_static_f = 0.5f * total_weight * (GP_LR * inv_wb);
    float fz_static_r = 0.5f * total_weight * (GP_LF * inv_wb);
    
    // Transferencia longitudinal por rueda (factor 0.5 para repartir el eje entre 2 ruedas)
    float dFz_lon_wheel = 0.5f * (GP_MASS * ax * GP_H_CG * inv_wb);

    // Transferencia lateral por rueda según el ancho de vía de cada eje
    float dFz_lat_f = 0.5f * (total_weight * (GP_LR * inv_wb) * (ay / GP_GRAVITY) * GP_H_CG) / GP_TRACK_F;
    float dFz_lat_r = 0.5f * (total_weight * (GP_LF * inv_wb) * (ay / GP_GRAVITY) * GP_H_CG) / GP_TRACK_R;
    
    // Carga aerodinámica (Downforce proporcional a v^2) en el eje trasero
    float downforce_rear_wheel = 0.25f * GP_AIR_DENSITY * (vx * vx) * GP_AERO_CL_REAR * GP_AERO_AREA;
    
    float fz_raw[4];
    fz_raw[GP_FL] = fz_static_f - dFz_lon_wheel + dFz_lat_f;
    fz_raw[GP_FR] = fz_static_f - dFz_lon_wheel - dFz_lat_f;
    fz_raw[GP_RL] = fz_static_r + dFz_lon_wheel + dFz_lat_r + downforce_rear_wheel;
    fz_raw[GP_RR] = fz_static_r + dFz_lon_wheel - dFz_lat_r + downforce_rear_wheel;
    
    for (int i = 0; i < 4; i++) {
        fz_out[i] = 50.0f + gp_softplus(fz_raw[i] - 50.0f);
    }
}

void gp_estimate_fy(float vx, float vy, float wz, float delta, const float fz[4], float fy_out[4]) {
    float vx_safe = fabsf(vx) + 0.5f;
    
    float alpha_f = delta - atan2f(vy + wz * GP_LF, vx_safe);
    float alpha_r =       - atan2f(vy - wz * GP_LR, vx_safe);
    
    float alpha_axle[2] = { alpha_f, alpha_r };
    float c_alpha[2]    = { GP_C_ALPHA_F, GP_C_ALPHA_R };
    int wheel_indices[2][2] = { {GP_FL, GP_FR}, {GP_RL, GP_RR} };

    for (int a = 0; a < 2; a++) {
        float alpha = alpha_axle[a];
        float sign_a = (alpha >= 0.0f) ? 1.0f : -1.0f;
        float abs_alpha = fabsf(alpha);

        for (int w = 0; w < 2; w++) {
            int idx = wheel_indices[a][w];
            float fz_w = fz[idx];

            // 1. Degresividad de fricción por carga vertical
            float mu_peak = GP_MU_NOM * (1.0f - GP_MU_DEG_FZ * (fz_w - 750.0f));
            mu_peak = GP_CLAMP(mu_peak, 0.8f, 1.8f);
            float mu_slide = GP_MU_SLIDE * (mu_peak / GP_MU_NOM);

            // 2. Parámetro de deformación del cepillo
            float theta = (0.5f * c_alpha[a] * abs_alpha) / (3.0f * mu_peak * fz_w + 1e-3f);

            float fy_w = 0.0f;
            if (theta < 1.0f) {
                // Zona elástica y de adherencia parcial
                fy_w = sign_a * (3.0f * mu_peak * fz_w * theta * (1.0f - theta + (theta * theta) / 3.0f));
            } else {
                // Zona de deslizamiento post-pico (decae suavemente hacia mu_slide * Fz)
                float decay = 1.0f / (1.0f + 0.5f * (theta - 1.0f));
                float mu_effective = mu_slide + (mu_peak - mu_slide) * decay;
                fy_w = sign_a * (mu_effective * fz_w);
            }

            fy_out[idx] = fy_w;
        }
    }
}

void gp_friction_ellipse_t_ub(const float fz[4], const float fy_est[4], float mu_est, float t_ub_out[4]) {
    float mu_safe = GP_CLAMP(mu_est, 0.4f, 2.0f);
    
    for (int i = 0; i < 4; i++) {
        float max_fy_capacity = mu_safe * fz[i];
        // Fx_sq_max = (mu*Fz)^2 - Fy^2
        float fx_sq_max = (max_fy_capacity * max_fy_capacity) - (fy_est[i] * fy_est[i]);
        
        // sqrt(softplus(x * 4) / 4)
        float fx_max = gp_softplus_sqrt(fx_sq_max * 4.0f) * 0.5f; 
        
        t_ub_out[i] = fx_max * GP_R_WHEEL;
    }
}

void gp_power_limited_t_ub(const float omega_wheel[4], float t_ub_out[4]) {
    for (int i = 0; i < 4; i++) {
        float omega_safe = gp_softplus(omega_wheel[i]);
        float t_power = GP_P_MAX_WHL / (omega_safe + 1e-3f);
        t_ub_out[i] = GP_CLAMP(t_power, 0.0f, GP_MAX_TRQ_ACT);
    }
}

void gp_power_limited_t_lb(const float omega_wheel[4], float p_max_charge_w, float t_lb_out[4]) {
    for (int i = 0; i < 4; i++) {
        float omega_safe = gp_softplus(fabsf(omega_wheel[i]));
        float t_power = p_max_charge_w / (omega_safe + 1e-3f);
        t_lb_out[i] = GP_CLAMP(t_power, 0.0f, GP_MAX_TRQ_ACT); 
    }
}

float gp_adaptive_k_us(const float fz[4]) {
    float fz_front_mean = 0.5f * (fz[GP_FL] + fz[GP_FR]);
    float fz_rear_mean  = 0.5f * (fz[GP_RL] + fz[GP_RR]);
    float delta_fz = fz_rear_mean - fz_front_mean;
    return GP_K_US * (1.0f + GP_K_US_FZ * delta_fz);
}

void gp_moment_arms(float arms_out[4]) {
    float hw_f = 0.5f * GP_TRACK_F / GP_R_WHEEL;
    float hw_r = 0.5f * GP_TRACK_R / GP_R_WHEEL;
    
    arms_out[GP_FL] = -hw_f;
    arms_out[GP_FR] =  hw_f;
    arms_out[GP_RL] = -hw_r;
    arms_out[GP_RR] =  hw_r;
}