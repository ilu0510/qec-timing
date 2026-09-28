# Paper Equations 

## Equation 1:
p_{idle} = 1 - \exp[-p \lambda \delta t] \approx p \lambda \delta t

## Equation 2:
p_{stab} = p
p_{data} = 1 - (1- p_{stab})(1- p_{idle}) \approx p (1 + \lambda \delta t)

## Equation 3:
p_{read} = b_{read} p

## Equation 4:
R = p_{L} / T

## Equation 5:
R = \frac{1}{\delta t} \frac{A}{d^{\beta}} (\frac{p}{p_{th}})^{(d+1)/2} (1+ \lambda \delta t)^{g (d+1)/2}

## Equation 6:
\Delta t^{*} = \operatorname*{argmin}_{\Delta t} R(\Delta t) = \frac{2}{\lambda (g(d+1)-2)} \approx \frac{2}{\lambda g d}

## Equation 9:
\lambda (t) = \lambda_{low} + (\lambda_{high} - \lambda_{low}) \chi (t), f = 1/T \int_0^T \chi (t) dt
\chi(t) \in \{0,1\} and r = \lambda_{high}/\lambda_{low}

## Equation A1:
D_{j, c} = s_{j, c} \oplus s_{j-1, c}

## Equation A2:
w_e = \log \frac{1-p_e}{p_e}

## Equation A3
p_L = 1/2 (1 - \exp(-2 \sum_j R(\delta t_j, t_j)\delta t_j))
t_j = \sum_{i<j} \Delta t_i

## Equation F3:
q_c(\lambda, \Delta t_j) = P(D_{j,c}=1 \mid \lambda, \Delta t_j)
= \frac{1}{2}\left[1-\prod_{e\in\partial(j,c)}(1-2p_e)\right]

## Equation F4:
\ell_j = \frac{1}{N_c}\sum_c
\left[
D_{j,c}\log\frac{q_c(\lambda_{\mathrm{high}},\Delta t_j)}
{q_c(\lambda_{\mathrm{low}},\Delta t_j)}
+
(1-D_{j,c})\log\frac{1-q_c(\lambda_{\mathrm{high}},\Delta t_j)}
{1-q_c(\lambda_{\mathrm{low}},\Delta t_j)}
\right]

## Equation F5:
\bar{\ell}^{(W)} = \frac{1}{W} \sum_{i=0}^{W-1} \ell_{j-i}.

## Equation F6:
\bar{\ell}^{(W)}



# Validation Checks

## Equation 7 / D9 (improvement of optimal over fixed interval)
\Gamma_{opt} = \frac{R(\Delta t)}{R(\Delta t^*)} = \frac{(1+\lambda\Delta t)^{\alpha}}{\lambda \Delta t K_\alpha},
\quad \alpha = g\frac{d+1}{2}, \quad K_\alpha = \frac{\alpha^\alpha}{(\alpha-1)^{\alpha-1}}

## Equation 10 (ideal adaptive advantage, large r)
\Gamma_{adapt} \approx \frac{r}{e \log r}, \quad f^*(r) \approx \frac{1}{r \log r}
(exact: \Gamma_{adapt}(10) \approx 1.86, \Gamma_{adapt}(20) \approx 2.73)

## Fitted constants
\beta = 2, g = 0.8 (b_read = 1, 2), g = 0.85 (b_read = 0.5)
p_th: 0.036 (b_read = 0.5), 0.029 (b_read = 1), 0.023 (b_read = 2)

## Figure parameters
Fig 2 / Fig 9: \lambda = 1, p = 0.015, b_read = 1, T = 200, d = 11..27; expect 1/\Delta t^* \approx 0.402 d + z
Fig 5a: d=5, p=0.005, T=100, A=0.75 | 5b: d=7, p=0.01, T=200, A=1.0 | 5c: d=11, p=0.015, T=500, A=1.4
Fig 5d: d=15, p=0.01, T=1000, A=1.5 | 5e: d=21, p=0.015, T=1000, A=1.9
Fig 3, d=7:  T=2000, p=0.015, \lambda_low=0.03, \lambda_high=2, \Delta t_low=1.7, \Delta t_high=0.25, f=0.064, W=7, \theta^*=0.045 -> \Gamma_{adapt} \approx 1.26
Fig 3, d=15: T=8000, p=0.01,  \lambda_low=0.03, \lambda_high=3, \Delta t_low=0.57, \Delta t_high=0.061, f=0.029, W=5, \theta^*=0.025 -> \Gamma_{adapt} \approx 1.95
Appendix C thresholds: \lambda = 0, 3d rounds