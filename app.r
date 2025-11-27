# Shiny app para los modelos SEIR determinista y estocástico
library(shiny)
library(shinythemes)
library(ggplot2)
library(plotly)
library(dplyr)
library(tidyr)
library(deSolve)
library(DT)
library(scales)

# -------------------------------
# Funciones de simulación basadas en los archivos originales
# -------------------------------
run_deterministic <- function(N, I0, E0, R0_init, gamma, sigma, niu, v, R0, t_max) {
  beta <- R0 * gamma
  S0 <- N - I0
  init <- c(S = S0, E = E0, I = I0, R = R0_init)
  t <- seq(0, t_max, 1)

  ed.sol <- function(t, state, parms) {
    with(as.list(state), {
      dxdt <- rep(0, length(state))
      dxdt[1] <- (niu * (N - S)) - (beta * I * S / N) - (v * S)
      dxdt[2] <- (beta * I * S / N) - ((niu + sigma) * E)
      dxdt[3] <- (sigma * E) - ((niu + gamma) * I)
      dxdt[4] <- gamma * I - (niu * R) + (v * S)
      list(dxdt)
    })
  }

  sol <- ode(y = init, times = t, func = ed.sol, parms = NULL)
  sol_df <- as.data.frame(sol) |> mutate(Modelo = "Determinista")
  sol_df
}

run_stochastic <- function(N, I0, E0, R0_init, gamma_base, sigma_base, niu, v, R0, t_max, dt,
                           periodo, pico_dia, fuerza_estacion, volatilidad_beta, volatilidad_bio,
                           n_simulaciones) {
  beta_base <- R0 * gamma_base

  sim_one <- function(id_sim) {
    t_vec <- seq(0, t_max, by = dt)
    n_steps <- length(t_vec)

    S <- numeric(n_steps); S[1] <- N - I0
    E <- numeric(n_steps); E[1] <- E0
    I <- numeric(n_steps); I[1] <- I0
    R <- numeric(n_steps); R[1] <- R0_init
    Beta_track <- numeric(n_steps); Beta_track[1] <- beta_base

    ruido_beta_actual <- rlnorm(1, 0, volatilidad_beta)
    ruido_bio_actual  <- rlnorm(1, 0, volatilidad_bio)

    for (i in 1:(n_steps - 1)) {
      S_now <- S[i]; E_now <- E[i]; I_now <- I[i]; R_now <- R[i]
      t_actual <- t_vec[i]

      if (floor(t_actual) > floor(t_actual - dt)) {
        ruido_beta_actual <- rlnorm(1, 0, volatilidad_beta)
        ruido_bio_actual  <- rlnorm(1, 0, volatilidad_bio)
      }

      factor_estacional <- 1 + fuerza_estacion * cos(2 * pi * (t_actual - pico_dia) / periodo)
      beta_now <- beta_base * factor_estacional * ruido_beta_actual
      Beta_track[i + 1] <- beta_now

      sigma_now <- sigma_base * ruido_bio_actual
      gamma_now <- gamma_base * ruido_bio_actual

      p_infec    <- 1 - exp(-(max(0, beta_now) * I_now / N) * dt)
      p_progress <- 1 - exp(-(max(0, sigma_now) * dt))
      p_recover  <- 1 - exp(-(max(0, gamma_now) * dt))
      p_death    <- 1 - exp(-niu * dt)
      p_vac      <- 1 - exp(-v * dt)

      new_SE <- rbinom(1, size = S_now, prob = p_infec)
      rem_S <- S_now - new_SE
      new_SR_vac <- rbinom(1, size = rem_S, prob = p_vac)
      rem_S <- rem_S - new_SR_vac
      death_S <- rbinom(1, size = rem_S, prob = p_death)

      new_EI   <- rbinom(1, size = E_now, prob = p_progress)
      death_E  <- rbinom(1, size = E_now - new_EI, prob = p_death)

      new_IR   <- rbinom(1, size = I_now, prob = p_recover)
      death_I  <- rbinom(1, size = I_now - new_IR, prob = p_death)

      death_R  <- rbinom(1, size = R_now, prob = p_death)
      new_births <- rpois(1, lambda = niu * N * dt)

      S[i + 1] <- max(0, S_now - new_SE - new_SR_vac - death_S + new_births)
      E[i + 1] <- max(0, E_now + new_SE - new_EI - death_E)
      I[i + 1] <- max(0, I_now + new_EI - new_IR - death_I)
      R[i + 1] <- max(0, R_now + new_IR + new_SR_vac - death_R)
    }

    data.frame(time = t_vec, S = S, E = E, I = I, R = R, Beta_Real = Beta_track,
               Sim = paste0("Sim ", id_sim), Modelo = "Estocástico")
  }

  sims <- lapply(seq_len(n_simulaciones), sim_one)
  bind_rows(sims)
}

# -------------------------------
# UI
# -------------------------------
header_css <- "
.navbar-default {background-color:#0c1b2a !important; border-color:#0c1b2a;}
.navbar-default .navbar-brand, .navbar-default .navbar-nav>li>a {color:#e5ecf4 !important; font-weight:600;}
.navbar-default .navbar-nav>li>a:hover {color:#64b5ff !important;}
body {background-color:#0b1624; color:#e5ecf4; padding-top:70px;}
.panel, .well {background-color:#0f2033; border-color:#10304d; color:#d9e3f3;}
.table {color:#d9e3f3;}
.btn-primary {background-color:#1e88e5; border-color:#1e88e5;}
.btn-primary:hover {background-color:#42a5f5;}
.param-card {border:1px solid #10304d; border-radius:6px; margin-bottom:15px;}
.param-card .panel-heading {cursor:pointer; display:flex; justify-content:space-between; align-items:center; background-color:#0f2a44; color:#e5ecf4; border-bottom:1px solid #10304d; padding:10px 12px;}
.param-card .panel-body {background-color:#0f2033;}
.param-card .close-card {color:#9fb3c8; font-size:16px; margin-left:10px;}
.param-card .close-card:hover {color:#ffffff;}
.simulate-box {text-align:center; margin:10px 0 20px;}
hr {border-top:1px solid #14406b;}
"

parameter_card <- function(card_id, title, body_content) {
  collapse_id <- paste0(card_id, "-body")
  div(
    class = "param-card panel panel-default", id = card_id,
    div(class = "panel-heading",
        span(title),
        span(
          tags$a(href = paste0("#", collapse_id), `data-toggle` = "collapse",
                 class = "text-info", icon("chevron-down")),
          tags$a(href = "#", class = "close-card pull-right", icon("xmark"))
        )
    ),
    div(id = collapse_id, class = "panel-collapse collapse in",
        div(class = "panel-body", body_content)
    )
  )
}

ui <- tagList(
  tags$head(
    tags$style(HTML(header_css)),
    tags$link(rel = "stylesheet",
              href = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css"),
    tags$script(HTML("$(document).on('click', '.close-card', function(e){e.preventDefault();$(this).closest('.param-card').remove();});"))
  ),
  navbarPage(
    title = "Modelos SEIR | Shiny",
    theme = shinytheme("cyborg"),
    windowTitle = "SEIR Dashboard",
    id = "top_nav",
    inverse = TRUE,

    tabPanel(
      "Modelos",
      fluidPage(
        br(),
        fluidRow(
          column(12,
                 tabsetPanel(
                   tabPanel("Curvas SEIR",
                           plotlyOutput("plot_series", height = "620px")),
                   tabPanel("Curva animada",
                            plotlyOutput("anim_plot", height = "620px")),
                   tabPanel("Tasa de contagio (β)",
                            plotlyOutput("beta_plot", height = "620px"))
                 )
          )
        ),
        hr(),
        fluidRow(
          column(4,
                 parameter_card(
                   "general-card",
                   tagList(icon("gears"), " Parámetros generales"),
                   tagList(
                     sliderInput("population", "Población total (N)",
                                 min = 1e5, max = 2e6, value = 1500000,
                                 step = 50000, sep = ",")
                   )
                 )
          ),
          column(4,
                 parameter_card(
                   "det-card",
                   tagList(icon("chart-line"), " Modelo determinista"),
                   tagList(
                     sliderInput("det_R0", "R0", min = 1, max = 6, value = 3.5, step = 0.1),
                     sliderInput("det_gamma", "Gamma (recuperación)",
                                 min = 0.05, max = 0.3, value = 0.143, step = 0.001),
                     sliderInput("det_sigma", "Sigma (latencia)",
                                 min = 0.05, max = 0.6, value = 0.25, step = 0.01),
                     sliderInput("det_niu", "Natalidad/Mortalidad (niu)",
                                 min = 0, max = 0.001, value = 0.0000245, step = 0.00001),
                     sliderInput("det_v", "Vacunación (v)",
                                 min = 0, max = 0.05, value = 0, step = 0.001),
                     sliderInput("det_I0", "Infectados iniciales (I0)",
                                 min = 1, max = 1000, value = 1),
                     sliderInput("det_E0", "Expuestos iniciales (E0)",
                                 min = 0, max = 2000, value = 0),
                     sliderInput("det_R0_init", "Recuperados iniciales (R0)",
                                 min = 0, max = 5000, value = 0),
                     sliderInput("det_tmax", "Días de simulación", min = 30, max = 365, value = 200, step = 5)
                   )
                 )
          ),
          column(4,
                 parameter_card(
                   "stoc-card",
                   tagList(icon("random"), " Modelo estocástico"),
                   tagList(
                     sliderInput("stoc_R0", "R0", min = 1, max = 6, value = 3.5, step = 0.1),
                     sliderInput("stoc_gamma", "Gamma base", min = 0.05, max = 0.3, value = 0.143, step = 0.001),
                     sliderInput("stoc_sigma", "Sigma base", min = 0.05, max = 0.6, value = 0.25, step = 0.01),
                     sliderInput("stoc_niu", "Natalidad/Mortalidad (niu)", min = 0, max = 0.001,
                                 value = 0.0000245, step = 0.00001),
                     sliderInput("stoc_v", "Vacunación (v)", min = 0, max = 0.05, value = 0, step = 0.001),
                     sliderInput("stoc_I0", "Infectados iniciales (I0)", min = 1, max = 2000, value = 5),
                     sliderInput("stoc_E0", "Expuestos iniciales (E0)", min = 0, max = 2000, value = 0),
                     sliderInput("stoc_R0_init", "Recuperados iniciales (R0)", min = 0, max = 5000, value = 0),
                     sliderInput("stoc_fuerza", "Fuerza estacional", min = 0, max = 1, value = 0.4, step = 0.05),
                     sliderInput("stoc_pico", "Pico estacional (día)", min = 0, max = 365, value = 75, step = 1),
                     sliderInput("stoc_periodo", "Periodo estacional (días)", min = 30, max = 365, value = 365, step = 5),
                     sliderInput("stoc_vol_beta", "Volatilidad contagio", min = 0, max = 1, value = 0.6, step = 0.05),
                     sliderInput("stoc_vol_bio", "Volatilidad recuperación", min = 0, max = 1, value = 0.3, step = 0.05),
                     sliderInput("stoc_dt", "Paso estocástico (dt)", min = 0.05, max = 1, value = 0.1, step = 0.05),
                     sliderInput("stoc_sims", "Número de simulaciones", min = 1, max = 30, value = 1, step = 1),
                     sliderInput("stoc_tmax", "Días de simulación", min = 30, max = 365, value = 200, step = 5)
                   )
                 )
          )
        ),
        fluidRow(
          column(12,
                 div(class = "simulate-box",
                     actionButton("simulate", label = "Simular", icon = icon("play"),
                                  class = "btn-primary btn-lg")
                 )
          )
        ),
        hr(),
        fluidRow(
          column(6,
                 h4("Valores deterministas"),
                 DTOutput("table_det")
          ),
          column(6,
                 h4("Valores estocásticos (media)"),
                 DTOutput("table_stoc")
          )
        )
      )
    )
  )
)

# -------------------------------
# Server
# -------------------------------
server <- function(input, output, session) {

  det_data <- eventReactive(input$simulate, {
    run_deterministic(
      N = input$population,
      I0 = input$det_I0,
      E0 = input$det_E0,
      R0_init = input$det_R0_init,
      gamma = input$det_gamma,
      sigma = input$det_sigma,
      niu = input$det_niu,
      v = input$det_v,
      R0 = input$det_R0,
      t_max = input$det_tmax
    )
  }, ignoreNULL = FALSE)

  stoc_data <- eventReactive(input$simulate, {
    run_stochastic(
      N = input$population,
      I0 = input$stoc_I0,
      E0 = input$stoc_E0,
      R0_init = input$stoc_R0_init,
      gamma_base = input$stoc_gamma,
      sigma_base = input$stoc_sigma,
      niu = input$stoc_niu,
      v = input$stoc_v,
      R0 = input$stoc_R0,
      t_max = input$stoc_tmax,
      dt = input$stoc_dt,
      periodo = input$stoc_periodo,
      pico_dia = input$stoc_pico,
      fuerza_estacion = input$stoc_fuerza,
      volatilidad_beta = input$stoc_vol_beta,
      volatilidad_bio = input$stoc_vol_bio,
      n_simulaciones = input$stoc_sims
    )
  }, ignoreNULL = FALSE)

  combined_series <- reactive({
    det <- det_data() |>
      pivot_longer(cols = c(S, E, I, R), names_to = "Compartimento", values_to = "Valor")
    st <- stoc_data() |>
      group_by(time, Modelo) |>
      summarise(across(c(S, E, I, R), mean), .groups = "drop") |>
      pivot_longer(cols = c(S, E, I, R), names_to = "Compartimento", values_to = "Valor")

    bind_rows(det, st)
  })

  animated_series <- reactive({
    df <- combined_series()
    tiempos <- sort(unique(df$time))
    bind_rows(lapply(tiempos, function(tiempo) {
      df |> filter(time <= tiempo) |> mutate(frame = tiempo)
    }))
  })

  output$plot_series <- renderPlotly({
    df <- combined_series()
    plot_ly(df, x = ~time, y = ~Valor, color = ~interaction(Modelo, Compartimento),
            colors = c(
              "Determinista.S" = "#42a5f5",
              "Determinista.E" = "#ff9800",
              "Determinista.I" = "#ef5350",
              "Determinista.R" = "#66bb6a",
              "Estocástico.S"  = "#1976d2",
              "Estocástico.E"  = "#ffa726",
              "Estocástico.I"  = "#e53935",
              "Estocástico.R"  = "#43a047"
            ),
            type = "scatter", mode = "lines",
            hovertemplate = paste0("<b>%{color}</b><br>Día %{x}<br>Valor: %{y:,.0f}<extra></extra>")) %>%
      layout(title = list(text = "Trayectorias SEIR"),
             xaxis = list(title = "Tiempo (días)"),
             yaxis = list(title = "Individuos", rangemode = "tozero"),
             plot_bgcolor = "#0f2033", paper_bgcolor = "#0f2033",
             font = list(color = "#e5ecf4"))
  })

  output$anim_plot <- renderPlotly({
    df <- animated_series()
    plot_ly(df, x = ~time, y = ~Valor,
            color = ~interaction(Modelo, Compartimento),
            colors = c(
              "Determinista.S" = "#42a5f5",
              "Determinista.E" = "#ff9800",
              "Determinista.I" = "#ef5350",
              "Determinista.R" = "#66bb6a",
              "Estocástico.S"  = "#1976d2",
              "Estocástico.E"  = "#ffa726",
              "Estocástico.I"  = "#e53935",
              "Estocástico.R"  = "#43a047"
            ),
            frame = ~frame, type = "scatter", mode = "lines",
            hovertemplate = paste0("<b>%{color}</b><br>Día %{x}<br>Valor: %{y:,.0f}<extra></extra>")) %>%
      layout(title = list(text = "Evolución animada SEIR"),
             xaxis = list(title = "Tiempo (días)"),
             yaxis = list(title = "Individuos", rangemode = "tozero"),
             plot_bgcolor = "#0f2033", paper_bgcolor = "#0f2033",
             font = list(color = "#e5ecf4"),
             updatemenus = list(list(
               type = "buttons",
               direction = "left",
               x = 0.1, y = 1.1, showactive = FALSE,
               buttons = list(
                 list(method = "animate", args = list(NULL, list(frame = list(duration = 100, redraw = TRUE), fromcurrent = TRUE)),
                      label = "Reproducir"),
                 list(method = "animate", args = list(NULL, list(mode = "immediate", frame = list(duration = 0, redraw = TRUE), transition = list(duration = 0))),
                      label = "Pausar")
               )
             )))
  })

  output$beta_plot <- renderPlotly({
    st <- stoc_data()
    beta_df <- st |>
      group_by(time) |>
      summarise(Beta_prom = mean(Beta_Real), .groups = "drop") |>
      mutate(Tendencia = input$stoc_R0 * input$stoc_gamma *
               (1 + input$stoc_fuerza *
                  cos(2 * pi * (time - input$stoc_pico) / input$stoc_periodo)))

    p <- ggplot(beta_df, aes(x = time)) +
      geom_line(aes(y = Beta_prom), color = "#ce93d8", size = 0.8) +
      geom_line(aes(y = Tendencia), linetype = "dashed",
                color = "#90caf9", size = 0.8) +
      labs(title = "Tasa de contagio (β)",
           subtitle = "Promedio estocástico vs. tendencia estacional",
           y = "β", x = "Tiempo") +
      theme_minimal(base_size = 13) +
      theme(plot.background = element_rect(fill = "#0f2033", colour = "#0f2033"),
            panel.background = element_rect(fill = "#0f2033", colour = NA),
            text = element_text(colour = "#e5ecf4"))

    ggplotly(p, tooltip = c("x", "y"))
  })

  output$table_det <- renderDT({
    det_data() |>
      select(time, S, E, I, R) |>
      datatable(options = list(pageLength = 6, scrollX = TRUE))
  })

  output$table_stoc <- renderDT({
    stoc_data() |>
      group_by(time) |>
      summarise(across(c(S, E, I, R), mean), .groups = "drop") |>
      datatable(options = list(pageLength = 6, scrollX = TRUE))
  })
}

shinyApp(ui, server)
