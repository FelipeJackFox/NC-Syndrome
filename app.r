# Shiny app combinando los modelos SEIR determinista y estocástico
library(shiny)
library(shinythemes)
library(ggplot2)
library(plotly)
library(dplyr)
library(tidyr)
library(deSolve)
library(DT)
library(scales)
library(rlang)

# -------------------------------
# Utilidades de simulación
# -------------------------------
run_deterministic <- function(population, init_I, init_E, init_R, gamma, sigma, nu, vaccination, R0, days, step) {
  beta <- R0 * gamma
  N <- population
  init <- c(S = N - init_I - init_E - init_R, E = init_E, I = init_I, R = init_R)
  times <- seq(0, days, by = step)
  
  ed.sol <- function(t, state, parms) {
    with(as.list(state), {
      dS <- (nu * (N - S)) - (beta * I * S / N) - (vaccination * S)
      dE <- (beta * I * S / N) - ((nu + sigma) * E)
      dI <- (sigma * E) - ((nu + gamma) * I)
      dR <- gamma * I - (nu * R) + (vaccination * S)
      list(c(dS, dE, dI, dR))
    })
  }
  
  sol <- ode(y = init, times = times, func = ed.sol, parms = NULL)
  sol_df <- as.data.frame(sol) |>
    mutate(Modelo = "Determinista")
  sol_df
}

run_stochastic <- function(population, init_I, init_E, init_R, gamma_base, sigma_base, nu, vaccination, R0,
                           days, dt, simulations, season_period, season_peak, season_strength,
                           vol_beta, vol_bio) {
  N <- population
  beta_base <- R0 * gamma_base
  t_vec <- seq(0, days, by = dt)
  n_steps <- length(t_vec)
  
  sim_one <- function(id_sim) {
    S <- numeric(n_steps); S[1] <- N - init_I - init_E - init_R
    E <- numeric(n_steps); E[1] <- init_E
    I <- numeric(n_steps); I[1] <- init_I
    R <- numeric(n_steps); R[1] <- init_R
    beta_tr <- numeric(n_steps); beta_tr[1] <- beta_base
    
    ruido_beta <- rlnorm(1, 0, vol_beta)
    ruido_bio  <- rlnorm(1, 0, vol_bio)
    
    for (i in 1:(n_steps - 1)) {
      S_now <- S[i]; E_now <- E[i]; I_now <- I[i]; R_now <- R[i]
      t_now <- t_vec[i]
      
      if (floor(t_now) > floor(t_now - dt)) {
        ruido_beta <- rlnorm(1, 0, vol_beta)
        ruido_bio  <- rlnorm(1, 0, vol_bio)
      }
      
      factor_estacional <- 1 + season_strength * cos(2 * pi * (t_now - season_peak) / season_period)
      beta_now <- beta_base * factor_estacional * ruido_beta
      beta_tr[i + 1] <- beta_now
      
      sigma_now <- sigma_base * ruido_bio
      gamma_now <- gamma_base * ruido_bio
      
      p_infec    <- 1 - exp(-(max(0, beta_now) * I_now / N) * dt)
      p_progress <- 1 - exp(-(max(0, sigma_now) * dt))
      p_recover  <- 1 - exp(-(max(0, gamma_now) * dt))
      p_death    <- 1 - exp(-nu * dt)
      p_vac      <- 1 - exp(-vaccination * dt)
      
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
      new_births <- rpois(1, lambda = nu * N * dt)
      
      S[i + 1] <- max(0, S_now - new_SE - new_SR_vac - death_S + new_births)
      E[i + 1] <- max(0, E_now + new_SE - new_EI - death_E)
      I[i + 1] <- max(0, I_now + new_EI - new_IR - death_I)
      R[i + 1] <- max(0, R_now + new_IR + new_SR_vac - death_R)
    }
    
    data.frame(time = t_vec, S = S, E = E, I = I, R = R,
               Beta = beta_tr, Sim = paste0("Sim ", id_sim))
  }
  
  sims <- lapply(seq_len(simulations), sim_one)
  bind_rows(sims, .id = "sim_id") |>
    mutate(Modelo = "Estocástico")
}

accumulate_by <- function(dat, var) {
  var <- rlang::as_string(ensym(var))
  lvls <- sort(unique(dat[[var]]))
  bind_rows(lapply(seq_along(lvls), function(i) {
    dat %>%
      filter(.data[[var]] <= lvls[i]) %>%
      mutate(frame = lvls[i])
  }))
}

# -------------------------------
# UI helpers
# -------------------------------
badge <- function(text, color = "#1e88e5") {
  tags$span(text, class = "badge-pill px-2 py-1",
            style = paste0("background-color:", color, "; color:white;"))
}

header_css <- "
.navbar-default {background-color:#0c1b2a !important; border-color:#0c1b2a;}
.navbar-default .navbar-brand, .navbar-default .navbar-nav>li>a {color:#e5ecf4 !important; font-weight:600;}
.navbar-default .navbar-nav>li>a:hover {color:#64b5ff !important;}
body {background-color:#0b1624; color:#e5ecf4; padding-top:70px;}
.panel, .well {background-color:#0f2033; border-color:#10304d; color:#d9e3f3;}
.table {color:#d9e3f3;}
.btn-primary {background-color:#1e88e5; border-color:#1e88e5;}
.btn-primary:hover {background-color:#42a5f5;}
hr {border-top:1px solid #14406b;}
.info-card {background:#0f2033; border:1px solid #174b7a; border-radius:12px; padding:18px; box-shadow:0 10px 24px rgba(0,0,0,0.35);} 
.mini-box {background:#0f2033; border:1px solid #174b7a; border-radius:12px; padding:14px 16px; margin-bottom:12px; box-shadow:0 4px 12px rgba(0,0,0,0.25);} 
.mini-box h5 {margin-top:0; color:#90caf9;}
.mini-box .value {font-size:26px; font-weight:700; color:#e5ecf4;}
"

context_cards <- function() {
  fluidRow(
    column(4,
           div(class = "info-card",
               h4(icon("bullseye"), "Propósito"),
               p("Aplicar modelado epidemiológico SEIR en sus variantes determinista y estocástica para informar la toma de decisiones en el brote descrito en el documento de contexto."),
               badge("Resumen ejecutivo")
           )
    ),
    column(4,
           div(class = "info-card",
               h4(icon("sliders"), "Controles interactivos"),
               p("Ajusta parámetros demográficos, tasas de transmisión y ruido estocástico en tiempo real. Observa cómo cambian las curvas en segundos, con animación inicial y actualizaciones reactivas."),
               badge("Laboratorio en vivo")
           )
    ),
    column(4,
           div(class = "info-card",
               h4(icon("palette"), "Diseño"),
               p("Tema oscuro con acentos azul eléctrico, gráfico animado con base \"ggplot\" y tablas estilizadas para lectura rápida."),
               badge("UI avanzada")
           )
    )
  )
}

# -------------------------------
# UI
# -------------------------------
ui <- tagList(
  tags$head(
    tags$style(HTML(header_css)),
    tags$link(rel = "stylesheet",
              href = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css")
  ),
  navbarPage(
    title = "Modelos SEIR | Shiny",
    theme = shinytheme("cyborg"),
    windowTitle = "SEIR Dashboard",
    id = "top_nav",
    inverse = TRUE,
    
    # ---- TAB 1: CONTEXTO ----
    tabPanel(
      "Contexto",
      fluidPage(
        br(),
        h2(icon("file-powerpoint"), "Situación problema"),
        p("El panel imita una diapositiva tipo presentación: a la izquierda se muestra el PDF fuente y a la derecha un resumen curado para discusión."),
        context_cards(),
        br(),
        fluidRow(
          column(7,
                 tags$div(style = "border:1px solid #174b7a; border-radius:12px; overflow:hidden;",
                          tags$iframe(
                            src = "Situación Problema.pdf",
                            width = "100%", height = "780px",
                            style = "border:none; background:#0f2033;"
                          )
                 )
          ),
          column(5,
                 div(class = "info-card",
                     h4(icon("notes-medical"), "Puntos clave"),
                     tags$ul(
                       tags$li("Identificar curvas epidémicas bajo diferentes supuestos de intervención."),
                       tags$li("Comparar sensibilidad de los parámetros de contagio y recuperación."),
                       tags$li("Contrastar predicciones deterministas vs. escenarios estocásticos con estacionalidad."),
                       tags$li("Disponer de tablas exportables para reportes ejecutivos."),
                       tags$li("Asegurar que cada curva respete la masa poblacional y muestre el error de balance.")
                     ),
                     hr(),
                     h5(icon("lightbulb"), "Sugerencias"),
                     p("Activa el refresco automático para observar trayectorias estocásticas dinámicas. Superpone ambos modelos para evaluar robustez."),
                     badge("Tip de facilitación", "#42a5f5")
                 )
          )
        )
      )
    ),
    
    # ---- TAB 2: MODELOS ----
    tabPanel(
      "Modelos",
      fluidPage(
        br(),
        # fila superior: gráficas + resumen
        fluidRow(
          column(9,
                 tabsetPanel(
                   tabPanel("Curvas y animación",
                            plotlyOutput("plot_series", height = "640px"),
                            br(),
                            plotOutput("error_plot", height = "210px")
                   ),
                   tabPanel("Tasa de contagio (β)",
                            plotlyOutput("beta_plot", height = "440px"))
                 )
          ),
          column(3,
                 wellPanel(
                   h4(icon("gauge-high"), "Resumen"),
                   br(),
                   uiOutput("peak_inf"),
                   uiOutput("peak_day"),
                   uiOutput("final_R")
                 )
          )
        ),
        hr(),
        # fila de parámetros
        fluidRow(
          column(4,
                 wellPanel(
                   h4(icon("gear"), "Parámetros generales"),
                   sliderInput("population", "Población total",
                               min = 1e5, max = 2e6, value = 1500000,
                               step = 50000, sep = ","),
                   sliderInput("days", "Días de simulación",
                               min = 30, max = 365, value = 200, step = 5),
                   sliderInput("step", "Paso de tiempo (días)",
                               min = 0.1, max = 5, value = 1, step = 0.1),
                   checkboxGroupInput("compartments", "Compartimentos a mostrar",
                                      choices = c("S", "E", "I", "R"),
                                      selected = c("S", "E", "I", "R")),
                   checkboxGroupInput("modelo_sel", "Modelo a visualizar",
                                      choices = c("Determinista", "Estocástico"),
                                      selected = c("Determinista", "Estocástico")),
                   checkboxInput("auto_refresh", "Refresco automático (estocástico)", TRUE),
                   actionButton("simulate", label = "Simular ahora",
                                icon = icon("play"),
                                class = "btn-primary btn-block")
                 )
          ),
          column(4,
                 wellPanel(
                   h4(icon("chart-line"), "Determinista"),
                   sliderInput("det_R0", "R0", min = 1, max = 6, value = 3.5, step = 0.1),
                   sliderInput("det_gamma", "Gamma (recuperación)",
                               min = 0.05, max = 0.3, value = 0.143, step = 0.001),
                   sliderInput("det_sigma", "Sigma (latencia)",
                               min = 0.05, max = 0.6, value = 0.25, step = 0.01),
                   sliderInput("det_nu", "Natalidad/Mortalidad (nu)",
                               min = 0, max = 0.001, value = 0.0000245, step = 0.00001),
                   sliderInput("det_v", "Vacunación",
                               min = 0, max = 0.05, value = 0, step = 0.001),
                   sliderInput("det_I0", "Infectados iniciales",
                               min = 1, max = 1000, value = 1),
                   sliderInput("det_E0", "Expuestos iniciales",
                               min = 0, max = 2000, value = 0),
                   sliderInput("det_R0_init", "Recuperados iniciales",
                               min = 0, max = 5000, value = 0)
                 )
          ),
          column(4,
                 wellPanel(
                   h4(icon("random"), "Estocástico"),
                   sliderInput("stoc_R0", "R0", min = 1, max = 6, value = 3.5, step = 0.1),
                   sliderInput("stoc_gamma", "Gamma base",
                               min = 0.05, max = 0.3, value = 0.143, step = 0.001),
                   sliderInput("stoc_sigma", "Sigma base",
                               min = 0.05, max = 0.6, value = 0.25, step = 0.01),
                   sliderInput("stoc_nu", "Natalidad/Mortalidad (nu)",
                               min = 0, max = 0.001, value = 0.0000245, step = 0.00001),
                   sliderInput("stoc_v", "Vacunación",
                               min = 0, max = 0.05, value = 0, step = 0.001),
                   sliderInput("stoc_I0", "Infectados iniciales",
                               min = 1, max = 2000, value = 5),
                   sliderInput("stoc_E0", "Expuestos iniciales",
                               min = 0, max = 2000, value = 0),
                   sliderInput("stoc_R0_init", "Recuperados iniciales",
                               min = 0, max = 5000, value = 0),
                   sliderInput("stoc_season_strength", "Fuerza estacional",
                               min = 0, max = 1, value = 0.4, step = 0.05),
                   sliderInput("stoc_season_peak", "Pico estacional (día)",
                               min = 0, max = 365, value = 75, step = 1),
                   sliderInput("stoc_period", "Periodo estacional (días)",
                               min = 30, max = 365, value = 365, step = 5),
                   sliderInput("stoc_vol_beta", "Volatilidad contagio",
                               min = 0, max = 1, value = 0.6, step = 0.05),
                   sliderInput("stoc_vol_bio", "Volatilidad recuperación",
                               min = 0, max = 1, value = 0.3, step = 0.05),
                   sliderInput("stoc_dt", "Paso estocástico (dt)",
                               min = 0.05, max = 1, value = 0.1, step = 0.05),
                   sliderInput("stoc_sims", "Número de simulaciones",
                               min = 1, max = 30, value = 5, step = 1)
                 )
          )
        ),
        hr(),
        # fila de tablas
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
  
  refreshTrigger <- reactiveVal(0)
  autoTimer <- reactiveTimer(2000)
  
  observeEvent(input$simulate, {
    refreshTrigger(isolate(refreshTrigger()) + 1)
  })
  
  observe({
    if (isTRUE(input$auto_refresh)) {
      autoTimer()
      isolate(refreshTrigger(refreshTrigger() + 1))
    }
  })
  
  det_data <- reactive({
    req("Determinista" %in% input$modelo_sel)
    run_deterministic(
      population = input$population,
      init_I = input$det_I0,
      init_E = input$det_E0,
      init_R = input$det_R0_init,
      gamma = input$det_gamma,
      sigma = input$det_sigma,
      nu = input$det_nu,
      vaccination = input$det_v,
      R0 = input$det_R0,
      days = input$days,
      step = input$step
    )
  })
  
  stoc_data <- reactive({
    req("Estocástico" %in% input$modelo_sel)
    refreshTrigger()
    run_stochastic(
      population = input$population,
      init_I = input$stoc_I0,
      init_E = input$stoc_E0,
      init_R = input$stoc_R0_init,
      gamma_base = input$stoc_gamma,
      sigma_base = input$stoc_sigma,
      nu = input$stoc_nu,
      vaccination = input$stoc_v,
      R0 = input$stoc_R0,
      days = input$days,
      dt = input$stoc_dt,
      simulations = input$stoc_sims,
      season_period = input$stoc_period,
      season_peak = input$stoc_season_peak,
      season_strength = input$stoc_season_strength,
      vol_beta = input$stoc_vol_beta,
      vol_bio = input$stoc_vol_bio
    )
  })
  
  combined_series <- reactive({
    data_list <- list()
    if ("Determinista" %in% input$modelo_sel) {
      data_list[["Determinista"]] <- det_data() |>
        pivot_longer(cols = c(S, E, I, R),
                     names_to = "Compartimento", values_to = "Valor") |>
        mutate(Modelo = "Determinista", frame = round(time))
    }
    if ("Estocástico" %in% input$modelo_sel) {
      st <- stoc_data()
      st_mean <- st |>
        group_by(time, Modelo) |>
        summarise(across(c(S, E, I, R), mean), .groups = "drop")
      data_list[["Estocástico"]] <- st_mean |>
        pivot_longer(cols = c(S, E, I, R),
                     names_to = "Compartimento", values_to = "Valor") |>
        mutate(Modelo = "Estocástico", frame = round(time))
    }
    if (length(data_list) == 0) return(data.frame())
    bind_rows(data_list)
  })
  
  output$plot_series <- renderPlotly({
    df <- combined_series()
    req(nrow(df) > 0)
    df <- df |> filter(Compartimento %in% input$compartments)
    
    palette <- c(
      "Determinista.S" = "#42a5f5",
      "Determinista.E" = "#ff9800",
      "Determinista.I" = "#ef5350",
      "Determinista.R" = "#66bb6a",
      "Estocástico.S"  = "#1976d2",
      "Estocástico.E"  = "#ffa726",
      "Estocástico.I"  = "#e53935",
      "Estocástico.R"  = "#43a047"
    )
    
    df_anim <- df |>
      mutate(curve = interaction(Modelo, Compartimento),
             curve = factor(curve, levels = names(palette))) |>
      accumulate_by(frame)
    
    plot_ly(
      df_anim,
      x = ~time,
      y = ~Valor,
      color = ~curve,
      colors = unname(palette),
      frame = ~frame,
      type = "scatter",
      mode = "lines",
      hovertemplate = paste0(
        "<b>%{color}</b><br>Día %{x}<br>Valor: %{y:,.0f}<extra></extra>"
      )
    ) %>%
      layout(
        title = list(text = "Trayectorias SEIR"),
        xaxis = list(title = "Tiempo (días)"),
        yaxis = list(title = "Individuos", rangemode = "tozero"),
        plot_bgcolor = "#0f2033",
        paper_bgcolor = "#0f2033",
        font = list(color = "#e5ecf4"),
        legend = list(orientation = "h", y = -0.18)
      ) %>%
      animation_opts(frame = 120, easing = "linear",
                     redraw = FALSE, transition = 0)
  })
  
  output$beta_plot <- renderPlotly({
    req("Estocástico" %in% input$modelo_sel)
    st <- stoc_data()
    beta_df <- st |>
      group_by(time) |>
      summarise(Beta_prom = mean(Beta), .groups = "drop") |>
      mutate(Tendencia = input$stoc_R0 * input$stoc_gamma *
               (1 + input$stoc_season_strength *
                  cos(2 * pi * (time - input$stoc_season_peak) /
                        input$stoc_period)))
    
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
    req("Determinista" %in% input$modelo_sel)
    det_data() |>
      select(time, S, E, I, R) |>
      datatable(options = list(pageLength = 6, scrollX = TRUE))
  })
  
  output$table_stoc <- renderDT({
    req("Estocástico" %in% input$modelo_sel)
    stoc_data() |>
      group_by(time) |>
      summarise(across(c(S, E, I, R), mean), .groups = "drop") |>
      datatable(options = list(pageLength = 6, scrollX = TRUE))
  })
  
  output$error_plot <- renderPlot({
    req("Determinista" %in% input$modelo_sel)
    det <- det_data()
    det <- det |> mutate(error_balance = (S + E + I + R) - input$population)
    ggplot(det, aes(x = time, y = error_balance)) +
      geom_area(fill = "#ef5350", alpha = 0.35) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "white") +
      labs(title = "Error de balance poblacional (Determinista)",
           x = "Tiempo", y = "S+E+I+R - N") +
      theme_minimal(base_size = 12) +
      theme(plot.background = element_rect(fill = "#0f2033", colour = "#0f2033"),
            panel.background = element_rect(fill = "#0f2033", colour = NA),
            text = element_text(colour = "#e5ecf4"))
  })
  
  observeEvent(combined_series(), {
    df <- combined_series()
    req(nrow(df) > 0)
    
    inf_curve <- df |>
      filter(Compartimento == "I") |>
      arrange(desc(Valor)) |>
      slice(1)
    
    rec_final <- df |>
      filter(Compartimento == "R", time == max(time)) |>
      summarise(valor = mean(Valor), .groups = "drop")
    
    render_mini_box <- function(title, value, accent = "#1e88e5") {
      div(class = "mini-box",
          h5(title),
          span(class = "value", value),
          div(style = paste0(
            "height:6px; background:linear-gradient(90deg,",
            accent, ", #0f2033); border-radius:10px; margin-top:8px;"
          ))
      )
    }
    
    output$peak_inf <- renderUI({
      render_mini_box("Pico de infectados",
                      scales::comma(round(inf_curve$Valor)),
                      "#1e88e5")
    })
    output$peak_day <- renderUI({
      render_mini_box("Momento del pico",
                      paste0(round(inf_curve$time), " días"),
                      "#26c6da")
    })
    output$final_R <- renderUI({
      render_mini_box("Recuperados finales",
                      scales::comma(round(rec_final$valor)),
                      "#43a047")
    })
  })
}

shinyApp(ui, server)
