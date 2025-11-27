library(ggplot2)
library(tidyr)
library(dplyr)

# Parámetros
N = 1500000
I0 = 5
E0 = 0 
S0 = N - I0
R0_init = 0

gamma_base = 0.143       
Rr = 3.5            
beta_base = Rr * gamma_base   
sigma_base = 0.25        
niu = 0.0000245     
v = 0                

# Configuración Estocástica
periodo <- 365
pico_dia <- 75         
fuerza_estacion <- 0.4 
volatilidad_beta <- 0.60  # Ruido contagio
volatilidad_bio <- 0.30   # Ruido recuperación

# Tiempo
t_max = 200
dt = 0.1 


# Simulación
correr_simulacion <- function(id_sim) {
  
  t_vec <- seq(0, t_max, by=dt)
  n_steps <- length(t_vec)
  
  S <- numeric(n_steps); S[1] <- S0
  E <- numeric(n_steps); E[1] <- E0
  I <- numeric(n_steps); I[1] <- I0
  R <- numeric(n_steps); R[1] <- R0_init
  Beta_track <- numeric(n_steps); Beta_track[1] <- beta_base
  
  # Ruidos del primer día
  ruido_beta_actual <- rlnorm(1, 0, volatilidad_beta)
  ruido_bio_actual  <- rlnorm(1, 0, volatilidad_bio)
  
  for(i in 1:(n_steps-1)){
    S_now <- S[i]; E_now <- E[i]; I_now <- I[i]; R_now <- R[i]
    t_actual <- t_vec[i]
    
    # Cuando cambia un día se genera un nuevo ruido
    if (floor(t_actual) > floor(t_actual - dt)) {
      ruido_beta_actual <- rlnorm(1, 0, volatilidad_beta)
      ruido_bio_actual  <- rlnorm(1, 0, volatilidad_bio)
    }
    
    # Estacionalidad y ruido contagio
    factor_estacional <- 1 + fuerza_estacion * cos(2 * pi * (t_actual - pico_dia) / periodo)
    beta_now <- beta_base * factor_estacional * ruido_beta_actual
    Beta_track[i+1] <- beta_now
    
    # Ruido recuperación
    sigma_now <- sigma_base * ruido_bio_actual
    gamma_now <- gamma_base * ruido_bio_actual 
    
    # Probabilidades
    p_infec    = 1 - exp(- (max(0, beta_now) * I_now / N) * dt)
    p_progress = 1 - exp(- (max(0, sigma_now) * dt))
    p_recover  = 1 - exp(- (max(0, gamma_now) * dt))
    p_death    = 1 - exp(- niu * dt)
    p_vac      = 1 - exp(- v * dt)
    
    # Binomiales
    new_SE = rbinom(1, size = S_now, prob = p_infec)
    rem_S = S_now - new_SE
    new_SR_vac = rbinom(1, size = rem_S, prob = p_vac)
    rem_S = rem_S - new_SR_vac
    death_S = rbinom(1, size = rem_S, prob = p_death)
    
    new_EI = rbinom(1, size = E_now, prob = p_progress)
    death_E = rbinom(1, size = E_now - new_EI, prob = p_death)
    
    new_IR = rbinom(1, size = I_now, prob = p_recover)
    death_I = rbinom(1, size = I_now - new_IR, prob = p_death)
    
    death_R = rbinom(1, size = R_now, prob = p_death)
    new_births = rpois(1, lambda = niu * N * dt)
    
    # Actualizaciones
    S[i+1] = max(0, S_now - new_SE - new_SR_vac - death_S + new_births)
    E[i+1] = max(0, E_now + new_SE - new_EI - death_E)
    I[i+1] = max(0, I_now + new_EI - new_IR - death_I)
    R[i+1] = max(0, R_now + new_IR + new_SR_vac - death_R)
  }
  
  return(data.frame(time = t_vec, S=S, E=E, I=I, R=R, Beta_Real=Beta_track, Sim = as.factor(id_sim)))
}


n_simulaciones = 1

lista_simulaciones <- lapply(1:n_simulaciones, correr_simulacion)
datos_totales <- do.call(rbind, lista_simulaciones)

datos_long <- pivot_longer(datos_totales, cols = c(S, E, I, R), 
                           names_to = "Compartimento", values_to = "Valor")

# Curvas SEIR
p1 <- ggplot(datos_long, aes(x = time, y = Valor, color = Compartimento, group = interaction(Sim, Compartimento))) +
  geom_line(alpha = 0.6, lwd = 0.4) + 
  labs(
    title = "Modelo SEIR Estocástico",
    subtitle = "Con estacionalidad y ruido diario",
    x = "Tiempo", y = "Individuos"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("S"="blue", "E"="orange", "I"="red", "R"="darkgreen"))

# Beta (Ritmo de Infección)
data_beta <- subset(datos_totales, Sim == 1)

# Tendencia ideal
data_beta$Tendencia <- beta_base * (1 + fuerza_estacion * cos(2 * pi * (data_beta$time - pico_dia) / periodo))

p2 <- ggplot(data_beta, aes(x = time)) +
  
  # Beta real (con ruido)
  geom_line(aes(y = Beta_Real), color = "purple", alpha = 0.6, lwd = 0.5) +
  
  # La tendencia estacional teórica
  geom_line(aes(y = Tendencia), color = "black", linetype = "dashed", lwd = 0.8) +
  labs(
    title = "Variación de la Tasa de Infección (Beta)",
    subtitle = "Morado: Beta Real con Ruido | Negro: Tendencia Estacional (Sin Ruido)",
    x = "Tiempo",
    y = "Valor de Beta"
  ) +
  theme_minimal()

print(p2)
print(p1)