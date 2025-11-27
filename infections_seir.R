
# Librerías de trabajo.
library(deSolve)
library(ggplot2)
library(tidyr) # Para pivot_longer


# Código de definición y resolución ED.
ed.sol = function(t, state, parms)
{
  with (as.list(state),
        {
          dxdt = rep(0, length(state))
          
          # dS/dt
          dxdt[1] = (niu*(N-S)) -(beta * I * S / N) - (v*S)
          
          # dE/dt
          dxdt[2] = (beta * I * S / N) - ((niu + sigma) * E)
          
          # dI/dt
          dxdt[3] = (sigma * E) - ((niu + gamma) * I)
          
          # dR/dt
          dxdt[4] = gamma * I - (niu * R) + (v*S)
          
          return(list(dxdt))
        })
}


# Condiciones Iniciales.

N = 1500000
I0 = 1
E0 = 0 
S0= N - I0
R0 = 0

Rr = 3.5

beta = Rr * gamma
gamma = 0.143  # Periodo infeccioso de 10 días.
sigma = 0.25  # Tasa de latencia.
niu = 0.0000245
v = 0  # No existe vacuna.

t = seq(0, 200, 1)

init = c(S = S0, E = E0, I = I0, R = R0)


miOutput = ode(y=init, times=t, func = ed.sol, parms = NULL)

miOutput <- as.data.frame(miOutput) 
head(miOutput)



datos_long <- pivot_longer(miOutput,
                           cols = -time,
                           names_to = "Compartimento",
                           values_to = "Valor")

ggplot(datos_long, aes(x = time, y = Valor, color = Compartimento)) +
  geom_line(lwd = 1.1) +
  labs(
    title = "Simulación de Modelo SEIR",
    subtitle = paste("N =", N, "| R0 =", R0),
    x = "Tiempo (días)",
    y = "Número de Individuos",
    color = "Compartimentos"
  ) +
  scale_y_continuous(labels = scales::comma) + 
  theme_minimal() +
  scale_color_manual(values = c(
    "S" = "blue",
    "E" = "orange",
    "I" = "red",
    "R" = "darkgreen"
  ))



colores <- c(S = "blue", E = "orange", I = "red", R = "darkgreen")

plot(miOutput$time, miOutput$S, type = "l",
     col = colores["S"],
     lwd = 2,
     ylim = c(0, N),
     xlab = "Tiempo (días)",
     ylab = "Número de Individuos",
     main = "Simulación de Modelo SEIR")

lines(miOutput$time, miOutput$E, type = "l", col = colores["E"], lwd = 2)
lines(miOutput$time, miOutput$I, type = "l", col = colores["I"], lwd = 2)
lines(miOutput$time, miOutput$R, type = "l", col = colores["R"], lwd = 2)

grid()

legend("right",
       legend = c("Susceptibles", "Expuestos", "Infectados", "Recuperados"),
       col = c(colores["S"], colores["E"], colores["I"], colores["R"]),
       lwd = 2,
       bty = "n")