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

    // 1. Ponderación objetivo de suavidad y regulación temporal
    const float t_bl_rl = (GP_W_REG * t_warmstart[GP_RL] + GP_W_SMOOTH * t_prev[GP_RL]) / h;
    const float t_bl_rr = (GP_W_REG * t_warmstart[GP_RR] + GP_W_SMOOTH * t_prev[GP_RR]) / h;

    const float lb_rl = t_lb[GP_RL], ub_rl = t_ub[GP_RL];
    const float lb_rr = t_lb[GP_RR], ub_rr = t_ub[GP_RR];

    // 2. Par total demandado por el piloto
    float t_target = fx_driver * GP_R_WHEEL;

    // Acotar la demanda al límite físico combinado de ambos actuadores
    const float t_total_min = lb_rl + lb_rr;
    const float t_total_max = ub_rl + ub_rr;
    t_target = GP_CLAMP(t_target, t_total_min, t_total_max);

    // 3. Solución libre proyectada sobre T_rl + T_rr = T_target
    float t_rl_sol = 0.5f * (t_target + t_bl_rl - t_bl_rr);
    float t_rr_sol = t_target - t_rl_sol;

    // 4. Árbol de decisión KKT exacto (garantiza Delta = 0.0 Nm si ambos saturan)
    if (t_rl_sol > ub_rl) {
        t_rl_sol = ub_rl;
        t_rr_sol = GP_CLAMP(t_target - t_rl_sol, lb_rr, ub_rr);
    } else if (t_rl_sol < lb_rl) {
        t_rl_sol = lb_rl;
        t_rr_sol = GP_CLAMP(t_target - t_rl_sol, lb_rr, ub_rr);
    } else if (t_rr_sol > ub_rr) {
        t_rr_sol = ub_rr;
        t_rl_sol = GP_CLAMP(t_target - t_rr_sol, lb_rl, ub_rl);
    } else if (t_rr_sol < lb_rr) {
        t_rr_sol = lb_rr;
        t_rl_sol = GP_CLAMP(t_target - t_rr_sol, lb_rl, ub_rl);
    }

    t_out[GP_FL] = 0.0f;
    t_out[GP_FR] = 0.0f;
    t_out[GP_RL] = t_rl_sol;
    t_out[GP_RR] = t_rr_sol;

    if (qp_residual != NULL) {
        *qp_residual = fabsf(a_eq * (t_rl_sol + t_rr_sol) - fx_driver);
    }
}