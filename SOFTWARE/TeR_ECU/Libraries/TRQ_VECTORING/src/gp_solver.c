#include "gp_solver.h"
#include <stddef.h>

void gp_nominal_allocation(float fx_driver, float mz_target, float t_nom_out[4]) {
    float arms[4];
    gp_moment_arms(arms);

    float t_fx = fx_driver * GP_R_WHEEL * 0.5f;
    float denom = (arms[GP_RL] * arms[GP_RL]) + (arms[GP_RR] * arms[GP_RR]);
    if (denom < 1e-6f) denom = 1e-6f;

    float t_mz_rl = (arms[GP_RL] * mz_target) / denom;
    float t_mz_rr = (arms[GP_RR] * mz_target) / denom;

    t_nom_out[GP_FL] = 0.0f;
    t_nom_out[GP_FR] = 0.0f;
    t_nom_out[GP_RL] = t_fx + t_mz_rl;
    t_nom_out[GP_RR] = t_fx + t_mz_rr;
}

void gp_qp_solve_rwd(
    const float t_warmstart[4],
    const float t_prev[4],
    float fx_driver,
    const float t_lb[4],
    const float t_ub[4],
    float alpha_qp,       
    float* lam_prev_ptr,  
    float t_out[4],
    float* qp_residual
) {
    float h = GP_W_REG + GP_W_SMOOTH;
    float a_eq = 1.0f / GP_R_WHEEL;
    float b_eq = fx_driver;

    float t_blend_rl = (GP_W_REG * t_warmstart[GP_RL] + GP_W_SMOOTH * t_prev[GP_RL]) / h;
    float t_blend_rr = (GP_W_REG * t_warmstart[GP_RR] + GP_W_SMOOTH * t_prev[GP_RR]) / h;

    float t_rl = GP_CLAMP(t_warmstart[GP_RL], t_lb[GP_RL], t_ub[GP_RL]);
    float t_rr = GP_CLAMP(t_warmstart[GP_RR], t_lb[GP_RR], t_ub[GP_RR]);
    
    float lam = *lam_prev_ptr;

    for (int i = 0; i < GP_QP_ITER; i++) {
        float viol = a_eq * (t_rl + t_rr) - b_eq;
        float g_rl = h * (t_rl - t_blend_rl) + a_eq * (lam + GP_RHO_AL * viol);
        float g_rr = h * (t_rr - t_blend_rr) + a_eq * (lam + GP_RHO_AL * viol);

        t_rl = GP_CLAMP(t_rl - alpha_qp * g_rl, t_lb[GP_RL], t_ub[GP_RL]);
        t_rr = GP_CLAMP(t_rr - alpha_qp * g_rr, t_lb[GP_RR], t_ub[GP_RR]);

        lam = lam + GP_RHO_AL * (a_eq * (t_rl + t_rr) - b_eq);
    }

    *lam_prev_ptr = GP_CLAMP(lam, -5000.0f, 5000.0f);

    t_out[GP_FL] = 0.0f;
    t_out[GP_FR] = 0.0f;
    t_out[GP_RL] = t_rl;
    t_out[GP_RR] = t_rr;

    if (qp_residual != NULL) {
        *qp_residual = fabsf(a_eq * (t_rl + t_rr) - b_eq);
    }
}

void gp_qp_solve_rwd_closedform(
    const float t_warmstart[4], const float t_prev[4], float fx_driver,
    const float t_lb[4], const float t_ub[4], float t_out[4], float* qp_residual
) {
    const float h    = GP_W_REG + GP_W_SMOOTH;
    const float a_eq = 1.0f / GP_R_WHEEL;
    const float b_eq = fx_driver;

  
    const float GP_SAT_SOFTNESS = 3.0f; // Nm

    const float t_bl_rl = (GP_W_REG * t_warmstart[GP_RL] + GP_W_SMOOTH * t_prev[GP_RL]) / h;
    const float t_bl_rr = (GP_W_REG * t_warmstart[GP_RR] + GP_W_SMOOTH * t_prev[GP_RR]) / h;

    const float lb_rl = t_lb[GP_RL], ub_rl = t_ub[GP_RL];
    const float lb_rr = t_lb[GP_RR], ub_rr = t_ub[GP_RR];

   
    const float lam = h * (a_eq * (t_bl_rl + t_bl_rr) - b_eq) / (2.0f * a_eq * a_eq);
    const float t_rl_free = t_bl_rl - lam * a_eq / h;
    const float t_rr_free = t_bl_rr - lam * a_eq / h;

   
    const float margin_rl = fminf(t_rl_free - lb_rl, ub_rl - t_rl_free);
    const float margin_rr = fminf(t_rr_free - lb_rr, ub_rr - t_rr_free);
    const float sat_rl = gp_sigmoid(-margin_rl / GP_SAT_SOFTNESS);
    const float sat_rr = gp_sigmoid(-margin_rr / GP_SAT_SOFTNESS);

    // --- Caso A: RL satura, RR resuelve.
    const float t_rl_A = GP_CLAMP(t_rl_free, lb_rl, ub_rl);
    const float t_rr_A = GP_CLAMP((b_eq / a_eq) - t_rl_A, lb_rr, ub_rr);

    // --- CAso B: RR satura, RL resuelve (simetrico a A)
    const float t_rr_B = GP_CLAMP(t_rr_free, lb_rr, ub_rr);
    const float t_rl_B = GP_CLAMP((b_eq / a_eq) - t_rr_B, lb_rl, ub_rl);

    // --- Candidato "ambos": clampeo independiente, imposible imposible
    const float t_rl_clamp_only = GP_CLAMP(t_rl_free, lb_rl, ub_rl);
    const float t_rr_clamp_only = GP_CLAMP(t_rr_free, lb_rr, ub_rr);

    // --- Mezcla ponderada de los 4 candidatos
    const float w_free = (1.0f - sat_rl) * (1.0f - sat_rr);
    const float w_A    = sat_rl * (1.0f - sat_rr);
    const float w_B    = (1.0f - sat_rl) * sat_rr;
    const float w_both = sat_rl * sat_rr;

    const float t_rl = w_free * t_rl_free + w_A * t_rl_A + w_B * t_rl_B + w_both * t_rl_clamp_only;
    const float t_rr = w_free * t_rr_free + w_A * t_rr_A + w_B * t_rr_B + w_both * t_rr_clamp_only;

    t_out[GP_FL] = 0.0f;
    t_out[GP_FR] = 0.0f;
    t_out[GP_RL] = t_rl;
    t_out[GP_RR] = t_rr;

    if (qp_residual != NULL) {
        *qp_residual = fabsf(a_eq * (t_rl + t_rr) - b_eq);
    }
}