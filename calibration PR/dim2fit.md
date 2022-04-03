Let $f(\boldsymbol x; \boldsymbol \theta) = \theta_3 x_2 + \theta_2 x_1 + \theta_1 $ as a linear combination of the input variables $\boldsymbol x$, where the input variable are $\boldsymbol x=[x_1, x_2]$, $x_1$ and $x_2$ are the input from each dimension; and $\boldsymbol \theta = [\theta_1, \theta_2, \theta_3]$ are corresponding parameters in the linear combination function $f$. The samples from the input space are denoted by $\boldsymbol \xi, \boldsymbol \zeta$.
$$
\nabla_\theta f (\boldsymbol x) = \frac{\partial }{\partial \boldsymbol \theta}f(\boldsymbol x; \boldsymbol \theta) = \left[ \begin{align*} 1 \\ x_1 \\ x_2 \end{align*} \right]
$$
The matrix constructed by gradients of randomly selected two samples $\boldsymbol \xi, \boldsymbol \zeta$
$$
\left (\nabla_\theta f(\boldsymbol \xi)\right) \left (\nabla_\theta f(\boldsymbol \zeta)\right)^T = \left[ \begin{matrix} 1 \\ \xi_1 \\ \xi_2 \end{matrix} \right] \left[ \begin{matrix} 1 & \zeta_1 & \zeta_2 \end{matrix} \right] = \left[ \begin{matrix} 1 & \zeta_1 & \zeta_2 \\ \xi_1 & \xi_1\zeta_1 & \xi_1\zeta_2 \\ \xi_2 & \xi_2\zeta_1 & \xi_2\zeta_2 \end{matrix} \right]
$$
The correlation kernel function of two samples $\boldsymbol \xi, \boldsymbol \zeta$ is:
$$
k_0(\boldsymbol \xi, \boldsymbol \zeta) = \exp(-\frac{(\xi_1-\zeta_1)^2}{\gamma_1})\exp(-\frac{(\xi_2-\zeta_2)^2}{\gamma_2})
$$
where $\gamma_1, \gamma_2$ are the parameters of kernel function in dimension $x_1$ and $x_2$.
$$
\boldsymbol H_\theta = \int_\boldsymbol \Xi \int_\boldsymbol \Zeta \left (\nabla_\theta f(\boldsymbol \xi)\right) \left (\nabla_\theta f(\boldsymbol \zeta)\right)^T k_0(\boldsymbol \xi, \boldsymbol \zeta) d \boldsymbol \xi d\boldsymbol \zeta
$$

$$
\nabla_\theta f (\boldsymbol x) = \frac{\partial }{\partial \boldsymbol \theta}f(\boldsymbol x; \boldsymbol \theta) = \left[ \begin{align*} 1 \\ x_1 \\ x_2 \end{align*} \right]
$$
The matrix constructed by gradients of randomly selected two samples $\boldsymbol \xi, \boldsymbol \zeta$
$$
\left (\nabla_\theta f(\boldsymbol \xi)\right) \left (\nabla_\theta f(\boldsymbol \zeta)\right)^T = \left[ \begin{matrix} 1 \\ \xi_1 \\ \xi_2 \end{matrix} \right] \left[ \begin{matrix} 1 & \zeta_1 & \zeta_2 \end{matrix} \right] = \left[ \begin{matrix} 1 & \zeta_1 & \zeta_2 \\ \xi_1 & \xi_1\zeta_1 & \xi_1\zeta_2 \\ \xi_2 & \xi_2\zeta_1 & \xi_2\zeta_2 \end{matrix} \right]
$$
The correlation kernel function of two samples $\boldsymbol \xi, \boldsymbol \zeta$ is:
$$
k_0(\boldsymbol \xi, \boldsymbol \zeta) = \exp(-\frac{(\xi_1-\zeta_1)^2}{\gamma_1})\exp(-\frac{(\xi_2-\zeta_2)^2}{\gamma_2})
$$
where $\gamma_1, \gamma_2$ are the parameters of kernel function in dimension $x_1$ and $x_2$.
$$
\boldsymbol H_\theta = \int_\boldsymbol \Xi \int_\boldsymbol \Zeta \left (\nabla_\theta f(\boldsymbol \xi)\right) \left (\nabla_\theta f(\boldsymbol \zeta)\right)^T g(\boldsymbol \xi, \boldsymbol \zeta) d \boldsymbol \xi d\boldsymbol \zeta
$$

$$
\nabla_\theta f (\boldsymbol x) = \frac{\partial }{\partial \boldsymbol \theta}f(\boldsymbol x; \boldsymbol \theta) = \left[ \begin{align*} 1 \\ x_1 \\ x_2 \end{align*} \right]
$$
The matrix constructed by gradients of randomly selected two samples $\boldsymbol \xi, \boldsymbol \zeta$
$$
\left (\nabla_\theta f(\boldsymbol \xi)\right) \left (\nabla_\theta f(\boldsymbol \zeta)\right)^T = \left[ \begin{matrix} 1 \\ \xi_1 \\ \xi_2 \end{matrix} \right] \left[ \begin{matrix} 1 & \zeta_1 & \zeta_2 \end{matrix} \right] = \left[ \begin{matrix} 1 & \zeta_1 & \zeta_2 \\ \xi_1 & \xi_1\zeta_1 & \xi_1\zeta_2 \\ \xi_2 & \xi_2\zeta_1 & \xi_2\zeta_2 \end{matrix} \right]
$$
The correlation kernel function of two samples $\boldsymbol \xi, \boldsymbol \zeta$ is:
$$
k_0(\boldsymbol \xi, \boldsymbol \zeta) = \exp(-\frac{(\xi_1-\zeta_1)^2}{\gamma_1})\exp(-\frac{(\xi_2-\zeta_2)^2}{\gamma_2})
$$
where $\gamma_1, \gamma_2$ are the parameters of kernel function in dimension $x_1$ and $x_2$.
$$
\boldsymbol H_\theta = \int_\boldsymbol \Xi \int_\boldsymbol \Zeta \left (\nabla_\theta f(\boldsymbol \xi)\right) \left (\nabla_\theta f(\boldsymbol \zeta)\right)^T k(\boldsymbol \xi, \boldsymbol \zeta) d \boldsymbol \xi d\boldsymbol \zeta
$$
