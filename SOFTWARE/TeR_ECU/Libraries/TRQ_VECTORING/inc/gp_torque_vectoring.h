#ifndef GP_TORQUE_VECTORING_H
#define GP_TORQUE_VECTORING_H

#include <stdint.h>
#include "gp_vehicle_model.h"
#include "gp_solver.h"           
#include "gp_traction_control.h" 
#include "gp_ekf.h"  


// Limites de control 
#define GP_TV_MAX_MZ                1500.0f
#define GP_TV_WZ_I_MAX              200.0f
#define GP_TV_RATE_LIMIT            3252.3f
#define GP_TV_EMA_ALPHA             0.2f
#define GP_MAX_BRAKE_PRESSURE_BAR   50.0f

// Ruiudo y filtrado 
#define GP_STEER_DEADZONE_RAD       0.0087f  
#define GP_YAW_DEADZONE_RADS        0.0175f  
#define GP_ACCEL_LPF_TAU            0.0200f  

#if defined(__arm__) || defined(__ARM_ARCH)
extern volatile uint32_t g_tv_exec_cycles;
extern volatile float g_tv_exec_us;
#endif

typedef struct {
    uint8_t enable;
    float   max_total_trq;
    float   max_charge_power_w;
} gp_regen_limits_t;

typedef struct {
    float wz_int;
    float delta_prev;
    float t_qp_prev[4];
    float t_out_prev[4];
    tc_state_t tc;
    gp_ekf_t ekf;
    float vy_est;
    float alpha_qp;
    float lam_prev;
    float mz_sat_ratio;
    float vy_gps_last;
    float vy_gps_age_ms;
    
    float ax_filt;
    float ay_filt;
    float t_ub_rl_filt;
    float t_ub_rr_filt;
    float t_lb_rl_filt;   
    float t_lb_rr_filt;   
    float qp_residual;
} tv_state_t;

void gp_tv_init(tv_state_t* state);

void gp_tv_step(
    float fx_driver, float delta, float vx, float vy, float wz, 
    float ay, float ax, const float omega[4], float brake_norm, 
    float temp_inv_rl, float temp_inv_rr, float vy_gps, uint8_t gps_valid,
    const gp_regen_limits_t* regen,
    float dt, tv_state_t* state, float t_cmd_out[4]
);

#endif // GP_TORQUE_VECTORING_H