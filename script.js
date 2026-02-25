const DENSITY = 80;

const sliders = {
    phi: document.getElementById('slider-phi'),
    mu: document.getElementById('slider-mu'),
    eta: document.getElementById('slider-eta'),
    epsilon: document.getElementById('slider-epsilon'),
    sigma: document.getElementById('slider-sigma'),
    beta: document.getElementById('slider-beta'),
    interest: document.getElementById('slider-interest'),
    zscale: document.getElementById('slider-zscale'),
    gamma: document.getElementById('slider-gamma')
};

const vals = {
    phi: document.getElementById('val-phi'),
    mu: document.getElementById('val-mu'),
    eta: document.getElementById('val-eta'),
    epsilon: document.getElementById('val-epsilon'),
    sigma: document.getElementById('val-sigma'),
    beta: document.getElementById('val-beta'),
    interest: document.getElementById('val-interest'),
    zscale: document.getElementById('val-zscale'),
    gamma: document.getElementById('val-gamma')
};

const slidersExt = {
    phi: document.getElementById('slider-phi-ext'),
    mu: document.getElementById('slider-mu-ext'),
    theta: document.getElementById('theta'),
    chi: document.getElementById('chi'),
    zeta: document.getElementById('zeta'),
    alpha: document.getElementById('alpha'),
    eta: document.getElementById('slider-eta-ext'),
    epsilon: document.getElementById('slider-epsilon-ext'),
    sigma: document.getElementById('slider-sigma-ext'),
    beta: document.getElementById('slider-beta-ext'),
    interest: document.getElementById('slider-interest-ext'),
    zscale: document.getElementById('slider-zscale-ext'),
    gamma: document.getElementById('slider-gamma-ext')
};

const valsExt = {
    phi: document.getElementById('val-phi-ext'),
    mu: document.getElementById('val-mu-ext'),
    theta: document.getElementById('theta-input'),
    chi: document.getElementById('chi-input'),
    zeta: document.getElementById('zeta-input'),
    alpha: document.getElementById('alpha-input'),
    eta: document.getElementById('val-eta-ext'),
    epsilon: document.getElementById('val-epsilon-ext'),
    sigma: document.getElementById('val-sigma-ext'),
    beta: document.getElementById('val-beta-ext'),
    interest: document.getElementById('val-interest-ext'),
    zscale: document.getElementById('val-zscale-ext'),
    gamma: document.getElementById('val-gamma-ext')
};

let currentMode = 'vary_epsilon';
let currentYMode = 'fix_phi';
let currentPolicyMode = 'active';
let currentModeExt = 'vary_epsilon';
let currentYModeExt = 'fix_phi';
let currentPolicyModeExt = 'active';

document.querySelectorAll('input[name="shock-type"]').forEach(radio => {
    radio.addEventListener('change', (e) => {
        currentMode = e.target.value;
        document.getElementById('group-eta').style.display = currentMode === 'vary_epsilon' ? 'flex' : 'none';
        document.getElementById('group-epsilon').style.display = currentMode === 'vary_eta' ? 'flex' : 'none';
        updatePlot();
    });
});

document.querySelectorAll('input[name="y-axis-mode"]').forEach(radio => {
    radio.addEventListener('change', (e) => {
        currentYMode = e.target.value;
        document.getElementById('group-phi').style.display = currentYMode === 'fix_phi' ? 'flex' : 'none';
        document.getElementById('group-mu').style.display = currentYMode === 'fix_mu' ? 'flex' : 'none';
        updatePlot();
    });
});

document.querySelectorAll('input[name="policy-mode"]').forEach(radio => {
    radio.addEventListener('change', (e) => {
        currentPolicyMode = e.target.value;
        document.getElementById('group-interest').style.display = currentPolicyMode === 'custom' ? 'flex' : 'none';
        updatePlot();
    });
});

document.querySelectorAll('input[name="shock-type-ext"]').forEach(radio => {
    radio.addEventListener('change', (e) => {
        currentModeExt = e.target.value;
        document.getElementById('group-eta-ext').style.display = currentModeExt === 'vary_epsilon' ? 'flex' : 'none';
        document.getElementById('group-epsilon-ext').style.display = currentModeExt === 'vary_eta' ? 'flex' : 'none';
        updatePlotExt();
    });
});

document.querySelectorAll('input[name="y-axis-mode-ext"]').forEach(radio => {
    radio.addEventListener('change', (e) => {
        currentYModeExt = e.target.value;
        document.getElementById('group-phi-ext').style.display = currentYModeExt === 'fix_phi' ? 'flex' : 'none';
        document.getElementById('group-mu-ext').style.display = currentYModeExt === 'fix_mu' ? 'flex' : 'none';
        updatePlotExt();
    });
});

document.querySelectorAll('input[name="policy-mode-ext"]').forEach(radio => {
    radio.addEventListener('change', (e) => {
        currentPolicyModeExt = e.target.value;
        document.getElementById('group-interest-ext').style.display = currentPolicyModeExt === 'custom' ? 'flex' : 'none';
        updatePlotExt();
    });
});

for (let key in vals) {
    if (vals[key] && vals[key].tagName === 'SPAN') {
        sliders[key].addEventListener('input', () => {
            vals[key].textContent = sliders[key].value;
            updatePlot();
        });
    }
}
for (let key in valsExt) {
    if (valsExt[key] && valsExt[key].tagName === 'SPAN') {
        slidersExt[key].addEventListener('input', () => {
            valsExt[key].textContent = slidersExt[key].value;
            updatePlotExt();
        });
    }
}
for (let key in valsExt) {
    if (valsExt[key] && valsExt[key].tagName === 'INPUT') {
        slidersExt[key].addEventListener('input', () => {
            valsExt[key].value = slidersExt[key].value;
            updatePlotExt();
        });
        valsExt[key].addEventListener('input', () => {
            let val = parseFloat(valsExt[key].value);
            if (!isNaN(val)) { slidersExt[key].value = val; updatePlotExt(); }
        });
    }
}

// ========== CORE MATH UTILS ==========
function bisection(f, a, b, tol = 1e-6, maxIter = 500) {
    let fa = f(a); let fb = f(b);
    if (fa * fb > 0) return null;
    let c = a;
    for (let i = 0; i < maxIter; i++) {
        c = (a + b) / 2; let fc = f(c);
        if (Math.abs(fc) < tol || (b - a) / 2 < tol) return c;
        if (fc * fa < 0) { b = c; fb = fc; }
        else { a = c; fa = fc; }
    }
    return c;
}

const CES = (ca, cb, ph, ep) => {
    if (ca <= 0 && cb <= 0) return 0;
    if (Math.abs(ep - 1.0) < 1e-5) return Math.pow(ca, ph) * Math.pow(cb, 1 - ph);
    let exponentInner = (ep - 1) / ep;
    if (exponentInner < -5) {
        let v_a = Math.log(ph) / ep + exponentInner * Math.log(ca);
        let v_b = Math.log(1 - ph) / ep + exponentInner * Math.log(cb);
        let max_v = Math.max(v_a, v_b);
        let sum_exp = Math.exp(v_a - max_v) + Math.exp(v_b - max_v);
        let log_res = max_v + Math.log(sum_exp);
        return Math.exp(log_res * (ep / (ep - 1)));
    }
    let inner = Math.pow(ph, 1 / ep) * Math.pow(ca, (ep - 1) / ep) + Math.pow(1 - ph, 1 / ep) * Math.pow(cb, (ep - 1) / ep);
    if (inner <= 0) return 0;
    let exponent = ep / (ep - 1);
    let result = Math.pow(inner, exponent);
    if (result === 0 || !isFinite(result)) { return Math.exp(exponent * Math.log(inner)); }
    return result;
}

function getCenter(arr_lower, arr_upper, deltas_arr) {
    let max_score = -1; let best_idx = 0;
    for (let i = 0; i < deltas_arr.length; i++) {
        let diff = arr_upper[i] - arr_lower[i];
        if (i < deltas_arr.length * 0.05 || i > deltas_arr.length * 0.95) continue;
        if (diff > max_score) { max_score = diff; best_idx = i; }
    }
    if (max_score === -1) best_idx = Math.floor(deltas_arr.length / 2);
    return { x: deltas_arr[best_idx], y: (arr_upper[best_idx] + arr_lower[best_idx]) / 2, diff: arr_upper[best_idx] - arr_lower[best_idx] };
}

// ========== 2D HEATMAP GRID BUILDER ==========
function get_state_for_Na(n_a, phi, is_ext, prm) {
    let retrained = (is_ext && prm.zeta) ? prm.zeta * Math.max(0, phi - n_a) : 0;
    let n_pot_a = phi - retrained;
    let n_pot_b = (1 - phi) + retrained;

    let y_a = prm.YA(n_a);
    let p_a = 1.0 / prm.MPLA(n_a);

    let p_b = p_a * Math.pow(((1 - phi) / Math.max(1e-5, phi)) * (y_a / Math.max(1e-5, n_pot_b)), 1 / prm.ep);
    let n_b = n_pot_b;

    if (p_b < 1.0) {
        p_b = 1.0;
        n_b = ((1 - phi) / Math.max(1e-5, phi)) * y_a * Math.pow(p_a, prm.ep);
        if (n_b > n_pot_b) n_b = n_pot_b;
    }

    let p_1 = Math.pow(phi * Math.pow(p_a, 1 - prm.ep) + (1 - phi) * Math.pow(p_b, 1 - prm.ep), 1 / (1 - prm.ep));
    let c_1 = CES(y_a, n_b, Math.max(1e-5, phi), prm.ep);

    let c_c;
    if (is_ext) {
        let profits = Math.max(0, p_a * y_a - n_a);
        let unemployed_income = prm.theta * (n_a / Math.max(1e-5, n_pot_a));
        let constrained_nominal_income = unemployed_income + prm.chi * (profits / Math.max(1e-5, n_pot_a));
        c_c = Math.max(0, constrained_nominal_income / p_1);
    } else {
        c_c = n_a / (Math.max(1e-5, phi) * p_1);
    }

    let c_u = prm.C_2 * Math.pow(Math.max(1e-8, prm.beta * (1 + prm.interest) * (p_1 / prm.P_2)), -prm.sigma);

    let denom = (is_ext ? n_pot_a : Math.max(1e-5, phi)) * Math.max(1e-6, c_u - c_c);
    let mu = (c_u - c_1) / denom;

    return { mu: mu, y: c_1, uA: Math.max(0, n_pot_a - n_a), uB: Math.max(0, n_pot_b - n_b) };
}

function solve_state(phi, target_mu, is_ext, prm, mu_g) {
    if (target_mu <= mu_g) {
        return get_state_for_Na(phi, phi, is_ext, prm);
    }

    let low = 1e-6, high = phi;
    let best_state = get_state_for_Na(phi, phi, is_ext, prm);
    for (let iter = 0; iter < 35; iter++) {
        let mid = (low + high) / 2;
        let state = get_state_for_Na(mid, phi, is_ext, prm);
        best_state = state;
        if (state.mu < target_mu) {
            high = mid;
        } else {
            low = mid;
        }
    }
    return best_state;
}

// ========== BASELINE SOLVER ==========
function get_mus_for_phi_delta(phi, delta) {
    let ep, et;
    if (currentMode === 'vary_epsilon') {
        et = parseFloat(sliders.eta.value);
        ep = et / (delta * et + 1.0);
    } else {
        ep = parseFloat(sliders.epsilon.value);
        let inv_et = 1.0 / ep - delta;
        et = 1.0 / inv_et;
    }
    if (et <= 1.001) et = 1.001;
    if (ep <= 1e-5) ep = 1e-5;

    let sigma = parseFloat(sliders.sigma.value);
    let beta = parseFloat(sliders.beta.value);

    let interest;
    if (currentPolicyMode === 'active') {
        interest = 0.0;
    } else if (currentPolicyMode === 'passive') {
        interest = 1.0 / beta - 1.0;
    } else {
        interest = parseFloat(sliders.interest.value);
    }

    let zscale = parseFloat(sliders.zscale.value);
    let gamma = parseFloat(sliders.gamma.value);

    let A = Math.pow(1 - gamma, -1 / (et - 1));
    let Z = zscale * phi * Math.pow((1 - gamma) / gamma, 1 / (et - 1));

    let exponentOuter = et / (et - 1);
    let exponentInner = (et - 1) / et;
    let Z_term = Math.pow(gamma, 1 / et) * Math.pow(Z, exponentInner);
    let N_coeff = Math.pow(1 - gamma, 1 / et);

    const YA = (N) => {
        if (N <= 0) return 0;
        let inner = Z_term + N_coeff * Math.pow(N, exponentInner);
        if (inner <= 0) return 0;
        return A * Math.pow(inner, exponentOuter);
    };

    const MPLA = (N) => {
        if (N <= 0) return Infinity;
        let Y = YA(N);
        return Math.pow(A, exponentInner) * N_coeff * Math.pow(Y / N, 1 / et);
    };

    const obj_p2 = (L2) => {
        let Y_a_2 = YA(L2);
        return Y_a_2 / (1 - L2) - (phi / (1 - phi)) * Math.pow(MPLA(L2), ep);
    };
    let L2 = bisection(obj_p2, 1e-6, phi);
    if (!L2) L2 = phi / 2;

    let Y_A2 = YA(L2);
    let Y_B2 = 1 - L2;
    let P_A2 = 1.0 / MPLA(L2);
    let P_B2 = 1.0;
    let P_2 = Math.pow(phi * Math.pow(P_A2, 1 - ep) + (1 - phi) * Math.pow(P_B2, 1 - ep), 1 / (1 - ep));
    let C_2 = CES(Y_A2, Y_B2, phi, ep);

    let Y_A1_gap = YA(phi);
    let P_A1_gap = 1.0 / MPLA(phi);
    let P_B1_gap = P_A1_gap * Math.pow(Y_A1_gap / phi, 1 / ep);
    let P_1_gap = Math.pow(phi * Math.pow(P_A1_gap, 1 - ep) + (1 - phi) * Math.pow(P_B1_gap, 1 - ep), 1 / (1 - ep));
    let C_gap = CES(Y_A1_gap, 1 - phi, phi, ep);

    let C_c_gap = 1.0 / P_1_gap;
    let C_u_gap = C_2 * Math.pow(Math.max(1e-8, beta * (1 + interest) * (P_1_gap / P_2)), -sigma);
    let denom_g = phi * Math.max(1e-6, C_u_gap - C_c_gap);
    let mu_g = (C_u_gap - C_gap) / denom_g;

    let N_A1_cont = bisection((N) => { return YA(N) - phi * Math.pow(MPLA(N), ep); }, 1e-6, phi);
    if (!N_A1_cont) N_A1_cont = phi - 1e-3;

    let P_A1_cont = 1.0 / MPLA(N_A1_cont);
    let P_1_cont = Math.pow(phi * Math.pow(P_A1_cont, 1 - ep) + (1 - phi), 1 / (1 - ep));
    let C_cont = Math.pow(P_1_cont, -ep);
    let C_c_cont = N_A1_cont / (phi * P_1_cont);
    let C_u_cont = C_2 * Math.pow(Math.max(1e-8, beta * (1 + interest) * (P_1_cont / P_2)), -sigma);
    let denom_c = phi * Math.max(1e-6, C_u_cont - C_c_cont);
    let mu_c = (C_u_cont - C_cont) / denom_c;

    let N_A1_rec = bisection((N) => {
        let p_a = 1.0 / MPLA(N);
        let p = Math.pow(phi * Math.pow(p_a, 1 - ep) + (1 - phi), 1 / (1 - ep));
        return YA(N) - phi * Math.pow(MPLA(N), ep) * Math.pow(p, ep);
    }, 1e-6, phi, 1e-8, 1000);
    if (!N_A1_rec) N_A1_rec = phi - 1e-3;

    let P_A1_rec = 1.0 / MPLA(N_A1_rec);
    let P_1_rec = Math.pow(phi * Math.pow(P_A1_rec, 1 - ep) + (1 - phi), 1 / (1 - ep));
    let C_c_rec = N_A1_rec / (phi * P_1_rec);
    let C_u_rec = C_2 * Math.pow(Math.max(1e-8, beta * (1 + interest) * (P_1_rec / P_2)), -sigma);
    let denom_r = phi * Math.max(1e-6, C_u_rec - C_c_rec);
    let mu_r = (C_u_rec - 1.0) / denom_r;

    let prm = { YA, MPLA, ep, beta, interest, sigma, C_2: C_2, P_2: P_2 };
    return { mu_g, mu_c, mu_r, prm: prm };
}

// ========== EXTENSIONS SOLVER ==========
function get_mus_for_phi_delta_ext(phi, delta) {
    let ep, et;
    if (currentModeExt === 'vary_epsilon') {
        et = parseFloat(slidersExt.eta.value);
        ep = et / (delta * et + 1.0);
    } else {
        ep = parseFloat(slidersExt.epsilon.value);
        let inv_et = 1.0 / ep - delta;
        et = 1.0 / inv_et;
    }

    if (et <= 1.001) et = 1.001;
    if (ep <= 1e-5) ep = 1e-5;

    let sigma = parseFloat(slidersExt.sigma.value);
    let beta = parseFloat(slidersExt.beta.value);

    let interest;
    if (currentPolicyModeExt === 'active') {
        interest = 0.0;
    } else if (currentPolicyModeExt === 'passive') {
        interest = 1.0 / beta - 1.0;
    } else {
        interest = parseFloat(slidersExt.interest.value);
    }

    let gamma = parseFloat(slidersExt.gamma.value);

    let theta = parseFloat(slidersExt.theta.value);
    let chi = parseFloat(slidersExt.chi.value);
    let zeta = parseFloat(slidersExt.zeta.value);
    let alpha = parseFloat(slidersExt.alpha.value);
    let zscale = parseFloat(slidersExt.zscale.value);

    let A = Math.pow(1 - gamma, -1 / (et - 1));
    let Base_Z = phi * Math.pow((1 - gamma) / gamma, 1 / (et - 1));

    let exponentOuter = et / (et - 1);
    let exponentInner = (et - 1) / et;
    let N_coeff = Math.pow(1 - gamma, 1 / et);

    const make_YA = (Z_val) => {
        let Z_term = Math.pow(gamma, 1 / et) * Math.pow(Z_val, exponentInner);
        return (N) => {
            if (N <= 0) return 0;
            let inner = Z_term + N_coeff * Math.pow(N, exponentInner);
            if (inner <= 0) return 0;
            return A * Math.pow(inner, exponentOuter);
        };
    };

    const make_MPLA = (YA_func) => {
        return (N) => {
            if (N <= 0) return Infinity;
            let Y = YA_func(N);
            return Math.pow(A, exponentInner) * N_coeff * Math.pow(Y / N, 1 / et);
        };
    };

    let Z_2 = zscale * Base_Z;
    const YA2 = make_YA(Z_2);
    const MPLA2 = make_MPLA(YA2);

    const obj_p2 = (L2) => {
        let Y_a_2 = YA2(L2);
        return Y_a_2 / (1 - L2) - (phi / (1 - phi)) * Math.pow(MPLA2(L2), ep);
    };
    let L2 = bisection(obj_p2, 1e-6, phi);
    if (!L2) L2 = phi / 2;

    let Y_A2 = YA2(L2);
    let Y_B2 = 1 - L2;
    let P_A2 = 1.0 / MPLA2(L2);
    let P_B2 = 1.0;
    let P_2 = Math.pow(phi * Math.pow(P_A2, 1 - ep) + (1 - phi) * Math.pow(P_B2, 1 - ep), 1 / (1 - ep));
    let C_2 = CES(Y_A2, Y_B2, phi, ep);

    let Z_1 = alpha * zscale * Base_Z;
    const YA = make_YA(Z_1);
    const MPLA = make_MPLA(YA);

    let Y_A1_gap = YA(phi);
    let P_A1_gap = 1.0 / MPLA(phi);
    let P_B1_gap = P_A1_gap * Math.pow(((1 - phi) / phi) * (Y_A1_gap / (1 - phi)), 1 / ep);
    let P_1_gap = Math.pow(phi * Math.pow(P_A1_gap, 1 - ep) + (1 - phi) * Math.pow(P_B1_gap, 1 - ep), 1 / (1 - ep));
    let C_gap = CES(Y_A1_gap, 1 - phi, phi, ep);

    const get_Cc = (N_A, Y_A, P_A, P_1, cur_n_pot_a) => {
        let labor_income = 1.0 * N_A;
        let profits = Math.max(0, P_A * Y_A - labor_income);
        let unemployed_income = theta * (N_A / Math.max(1e-5, cur_n_pot_a));
        let constrained_nominal_income = unemployed_income + chi * (profits / Math.max(1e-5, cur_n_pot_a));
        return Math.max(0, constrained_nominal_income / P_1);
    };

    let C_c_gap = get_Cc(phi, Y_A1_gap, P_A1_gap, P_1_gap, phi);
    let C_u_gap = C_2 * Math.pow(Math.max(1e-8, beta * (1 + interest) * (P_1_gap / P_2)), -sigma);
    let denom_g = phi * Math.max(1e-6, C_u_gap - C_c_gap);
    let mu_g = (C_u_gap - C_gap) / denom_g;

    let N_A1_cont = bisection((N) => {
        let retrained = zeta * (phi - N);
        let n_sb = (1 - phi) + retrained;
        return YA(N) - (phi / (1 - phi)) * n_sb * Math.pow(MPLA(N), ep);
    }, 1e-6, phi);
    if (!N_A1_cont) N_A1_cont = phi - 1e-3;

    let retrained_cont = zeta * (phi - N_A1_cont);
    let n_pot_a_cont = phi - retrained_cont;
    let n_pot_b_cont = (1 - phi) + retrained_cont;

    let P_A1_cont = 1.0 / MPLA(N_A1_cont);
    let P_1_cont = Math.pow(phi * Math.pow(P_A1_cont, 1 - ep) + (1 - phi), 1 / (1 - ep));
    let C_cont = CES(YA(N_A1_cont), n_pot_b_cont, phi, ep);
    let C_c_cont = get_Cc(N_A1_cont, YA(N_A1_cont), P_A1_cont, P_1_cont, n_pot_a_cont);
    let C_u_cont = C_2 * Math.pow(Math.max(1e-8, beta * (1 + interest) * (P_1_cont / P_2)), -sigma);
    let denom_c = n_pot_a_cont * Math.max(1e-6, C_u_cont - C_c_cont);
    let mu_c = (C_u_cont - C_cont) / denom_c;

    let N_A1_rec = bisection((N) => {
        let p_a = 1.0 / MPLA(N);
        let p = Math.pow(phi * Math.pow(p_a, 1 - ep) + (1 - phi), 1 / (1 - ep));
        return YA(N) - phi * Math.pow(MPLA(N), ep) * Math.pow(p, ep);
    }, 1e-6, phi, 1e-8, 1000);
    if (!N_A1_rec) N_A1_rec = phi - 1e-3;

    let retrained_rec = zeta * (phi - N_A1_rec);
    let n_pot_a_rec = phi - retrained_rec;

    let P_A1_rec = 1.0 / MPLA(N_A1_rec);
    let P_1_rec = Math.pow(phi * Math.pow(P_A1_rec, 1 - ep) + (1 - phi), 1 / (1 - ep));
    let C_c_rec = get_Cc(N_A1_rec, YA(N_A1_rec), P_A1_rec, P_1_rec, n_pot_a_rec);
    let C_u_rec = C_2 * Math.pow(Math.max(1e-8, beta * (1 + interest) * (P_1_rec / P_2)), -sigma);
    let denom_r = n_pot_a_rec * Math.max(1e-6, C_u_rec - C_c_rec);
    let mu_r = (C_u_rec - 1.0) / denom_r;

    let prm = { YA, MPLA, ep, beta, interest, sigma, C_2: C_2, P_2: P_2, zeta: zeta, theta: theta, chi: chi };
    return { mu_g, mu_c, mu_r, prm: prm };
}

// ========== PLOT GENERATORS ==========
function drawPlotHelper(divId, getMusFunc, sliderMap, mode, yMode) {
    let deltas = [];
    let varphi_gap = []; let varphi_cont = []; let varphi_rec = [];
    let params_list = [];

    let log_delta_min = -2;
    if (mode === 'vary_epsilon') {
        let log_delta_max = 2.0;
        let step = (log_delta_max - log_delta_min) / DENSITY;
        for (let i = 0; i <= DENSITY; i++) { deltas.push(Math.pow(10, log_delta_min + i * step)); }
    } else {
        let ep = parseFloat(sliderMap.epsilon.value);
        let max_delta = 1.0 / ep;
        let delta_max = max_delta * 0.999;
        let step = (delta_max - 0.01) / DENSITY;
        for (let i = 0; i <= DENSITY; i++) { deltas.push(0.01 + i * step); }
    }

    let y_max = 1.0;
    for (let i = 0; i < deltas.length; i++) {
        let delta = deltas[i];
        let fg, fc, fr, prm;

        if (yMode === 'fix_phi') {
            let phi = parseFloat(sliderMap.phi.value);
            let res = getMusFunc(phi, delta);
            fg = phi * Math.max(0, Math.min(1, res.mu_g));
            fc = phi * Math.max(0, Math.min(1, res.mu_c));
            fr = phi * Math.max(0, Math.min(1, res.mu_r));
            prm = res.prm;
            y_max = phi;
        } else {
            let mu = parseFloat(sliderMap.mu.value);
            let solve_phi = (prop) => {
                let f = (phi) => getMusFunc(phi, delta)[prop] - mu;
                let left = 1e-4, right = 0.999;
                let fl = f(left), fr = f(right);
                if (fl * fr <= 0) return bisection(f, left, right);
                return fl > 0 ? 1.0 : 0.0;
            };
            let ans_g = solve_phi('mu_g');
            let ans_c = solve_phi('mu_c');
            let ans_r = solve_phi('mu_r');
            fg = mu * (ans_g !== null ? ans_g : 1.0);
            fc = mu * (ans_c !== null ? ans_c : 1.0);
            fr = mu * (ans_r !== null ? ans_r : 1.0);
            // Get params for the median phi roughly just for grid mapping (simplification)
            let median_phi = 0.5;
            prm = getMusFunc(median_phi, delta).prm;
            y_max = mu;
        }

        let final_g = Math.max(0.0, Math.min(y_max, fg));
        let final_c = Math.max(final_g, Math.min(y_max, fc));
        let final_r = Math.max(final_c, Math.min(y_max, fr));

        varphi_gap.push(final_g); varphi_cont.push(final_c); varphi_rec.push(final_r);
    }

    // Heatmap traces
    let num_y_points = 100;
    let varphi_grid = [];
    for (let j = 0; j <= num_y_points; j++) {
        varphi_grid.push((j / num_y_points) * y_max);
    }

    let Z_out = [];
    let Z_uA = [];
    let Z_uB = [];
    for (let j = 0; j <= num_y_points; j++) {
        Z_out.push(new Array(deltas.length).fill(0));
        Z_uA.push(new Array(deltas.length).fill(0));
        Z_uB.push(new Array(deltas.length).fill(0));
    }

    let is_ext = (divId === 'plot-div-ext');

    for (let i = 0; i < deltas.length; i++) {
        let delta = deltas[i];
        let prm_cache = null;
        let mu_g_cache = null;

        if (yMode === 'fix_phi') {
            let phi = parseFloat(sliderMap.phi.value);
            let res = getMusFunc(phi, delta);
            prm_cache = res.prm;
            mu_g_cache = res.mu_g;
        }

        for (let j = 0; j <= num_y_points; j++) {
            let varphi = varphi_grid[j];
            let phi, target_mu, prm, mu_g;

            if (yMode === 'fix_phi') {
                phi = parseFloat(sliderMap.phi.value);
                target_mu = varphi / Math.max(1e-4, phi);
                prm = prm_cache;
                mu_g = mu_g_cache;
            } else {
                target_mu = parseFloat(sliderMap.mu.value);
                phi = Math.max(1e-4, Math.min(0.999, varphi / Math.max(1e-4, target_mu)));
                let res = getMusFunc(phi, delta);
                prm = res.prm;
                mu_g = res.mu_g;
            }

            let state = solve_state(phi, target_mu, is_ext, prm, mu_g);
            Z_out[j][i] = state.y;
            Z_uA[j][i] = state.uA;
            Z_uB[j][i] = state.uB;
        }
    }

    let trace_bottom = { x: deltas, y: deltas.map(() => 0), mode: 'lines', line: { width: 0 }, showlegend: false, hoverinfo: 'skip' };
    let trace_gap_line = { x: deltas, y: varphi_gap, mode: 'lines', name: 'Output Gap', line: { color: 'black', width: 4, dash: 'solid' }, fill: 'none' };
    let trace_cont_line = { x: deltas, y: varphi_cont, mode: 'lines', name: 'Contagion', line: { color: 'black', width: 4, dash: 'dash' }, fill: 'none' };
    let trace_rec_line = { x: deltas, y: varphi_rec, mode: 'lines', name: 'Recession', line: { color: 'black', width: 4, dash: 'dot' }, fill: 'none' };
    let trace_top = { x: deltas, y: deltas.map(() => y_max), mode: 'lines', showlegend: false, line: { width: 0 } };

    let arr_bottom = deltas.map(() => 0);
    let arr_top = deltas.map(() => y_max);
    let c_fe = getCenter(arr_bottom, varphi_gap, deltas);
    let c_tu = getCenter(varphi_gap, varphi_cont, deltas);
    let c_ku = getCenter(varphi_cont, varphi_rec, deltas);
    let c_tr = getCenter(varphi_rec, arr_top, deltas);

    let xFn = (x) => mode === 'vary_epsilon' ? Math.log10(x) : x;
    let computed_annotations = [];
    if (c_fe.diff > 0.02) computed_annotations.push({ x: xFn(c_fe.x), y: c_fe.y, text: 'Full<br>Employment', showarrow: false, xanchor: 'left', font: { color: '#ffffff', size: 14 } });
    if (c_tu.diff > 0.02) computed_annotations.push({ x: xFn(c_tu.x), y: c_tu.y, text: 'Tech<br>Unemp.', showarrow: false, font: { color: '#ffffff', size: 14 } });
    if (c_ku.diff > 0.02) computed_annotations.push({ x: xFn(c_ku.x), y: c_ku.y, text: 'Keynesian<br>Unemp.', showarrow: false, font: { color: '#ffffff', size: 14 } });
    if (c_tr.diff > 0.02) computed_annotations.push({ x: xFn(c_tr.x), y: c_tr.y, text: 'Tech<br>Recession', showarrow: false, xanchor: 'right', font: { color: '#ffffff', size: 16, weight: 'bold' } });

    let xaxis_conf = { title: 'Sensitivity of Labor Displacement, δ = 1/ε - 1/η', color: '#111827', gridcolor: '#e5e7eb', zerolinecolor: '#e5e7eb', range: mode === 'vary_epsilon' ? [-2, 2.0] : [0, 1.0 / parseFloat(sliderMap.epsilon.value)] };
    if (mode === 'vary_epsilon') xaxis_conf.type = 'log';

    let layout_out = {
        title: { text: (is_ext ? 'Extensions Output' : 'Baseline Output'), font: { color: '#111827', size: 20 } },
        paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
        height: 500,
        xaxis: xaxis_conf,
        yaxis: { title: 'Economic Exposure, φ = ϕμ', color: '#111827', gridcolor: '#e5e7eb', zerolinecolor: '#111827', range: [0, y_max] },
        legend: { orientation: 'h', x: -0.05, y: -0.25, xanchor: 'left', yanchor: 'top', itemwidth: 30, bgcolor: 'rgba(255, 255, 255, 0.8)', bordercolor: '#e5e7eb', borderwidth: 1, font: { color: '#111827', size: 11 } },
        margin: { l: 60, r: 120, t: 80, b: 70 },
        annotations: computed_annotations
    };

    let layout_u = {
        title: { text: (is_ext ? 'Extensions Unemployment' : 'Baseline Unemployment'), font: { color: '#111827', size: 20 } },
        paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
        height: 500,
        xaxis: xaxis_conf,
        yaxis: { title: '', color: '#111827', gridcolor: '#e5e7eb', zerolinecolor: '#111827', range: [0, y_max] },
        showlegend: true,
        legend: { orientation: 'h', x: -0.05, y: -0.25, xanchor: 'left', yanchor: 'top', itemwidth: 20, bgcolor: 'rgba(255, 255, 255, 0.8)', bordercolor: '#e5e7eb', borderwidth: 1, font: { color: '#111827', size: 10, family: 'Inter, sans-serif' } },
        margin: { l: 60, r: 160, t: 80, b: 70 },
        annotations: computed_annotations
    };

    let trace_out_heatmap = {
        x: deltas,
        y: varphi_grid,
        z: Z_out,
        type: 'contour',
        colorscale: 'Portland',
        reversescale: true,
        showscale: true,
        zmin: 0.90,
        zmax: 1.25,
        hoverinfo: 'x+y+z',
        contours: {
            coloring: 'heatmap',
            showlines: false
        },
        colorbar: {
            title: 'Output',
            titleside: 'top'
        }
    };

    let trace_uA = {
        x: deltas,
        y: varphi_grid,
        z: Z_uA,
        type: 'contour',
        colorscale: [[0, 'rgba(255, 0, 0, 0)'], [0.1, 'rgba(255, 0, 0, 0.45)'], [1, 'rgba(255, 0, 0, 0.95)']],
        showscale: true,
        zmin: 0.0,
        zmax: 0.35,
        hoverinfo: 'x+y+z',
        contours: {
            coloring: 'heatmap',
            showlines: false
        },
        colorbar: {
            title: 'Sector A<br>Unemp.',
            titleside: 'top',
            x: 1.05
        }
    };

    let trace_uB = {
        x: deltas,
        y: varphi_grid,
        z: Z_uB,
        type: 'contour',
        colorscale: [[0, 'rgba(0, 0, 255, 0)'], [0.1, 'rgba(0, 0, 255, 0.45)'], [1, 'rgba(0, 0, 255, 0.95)']],
        showscale: true,
        zmin: 0.0,
        zmax: 0.35,
        hoverinfo: 'x+y+z',
        contours: {
            coloring: 'heatmap',
            showlines: false
        },
        colorbar: {
            title: 'Sector B<br>Unemp.',
            titleside: 'top',
            x: 1.30
        }
    };

    let lines = [trace_gap_line, trace_cont_line, trace_rec_line];

    Plotly.react(divId + '-output', [trace_out_heatmap, ...lines], layout_out, { responsive: true });
    Plotly.react(divId + '-unemp', [trace_uA, trace_uB, ...lines], layout_u, { responsive: true });
}

function updatePlot() {
    drawPlotHelper('plot-div', get_mus_for_phi_delta, sliders, currentMode, currentYMode);
}

function updatePlotExt() {
    drawPlotHelper('plot-div-ext', get_mus_for_phi_delta_ext, slidersExt, currentModeExt, currentYModeExt);
}

updatePlot();
updatePlotExt();
