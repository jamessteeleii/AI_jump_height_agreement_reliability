##### Function to prepare data
prepare_data <- function(file) {
  dat <- read_csv(file) |>
    janitor::clean_names() |>
    pivot_longer(5:22,
                 names_to = "method",
                 values_to = "jump_height") |>
    separate(method, c("method","jump", "trial", "session")) |>
    select(-jump) |>
    filter(method != "FDFT") |>
    mutate(trial = as.numeric(trial)) |>
    mutate(session = as.numeric(
      case_when(
        session == "WK1" ~ 1,
        session == "WK2" ~ 2
      )
    ))
}

##### Functions to fit agreement and reliability models
fit_agree_models <- function(data) {
  
  agree_models <- tibble(method = as.character(),
                         mean_bias = as.numeric(),
                         lower_ci_mean = as.numeric(), 
                         upper_ci_mean = as.numeric(), 
                         lower_ci_l_loa = as.numeric(), 
                         upper_ci_l_loa = as.numeric(),
                         lower_ci_u_loa = as.numeric(),
                         upper_ci_u_loa = as.numeric(),
                         l_loa = as.numeric(),
                         u_loa = as.numeric())
  
  set.seed(1988)
  
  for(i in c("jmc", "mj", "opto")) {
    
    data_wide <- data |>
      filter(method == "fdi" | method == i) |>
      mutate(method = case_when(
        method == "fdi" ~ "M1",
        method == "jmc" ~ "M2",
        method == "mj" ~ "M2",
        method == "opto" ~ "M2"
      )) |>
      pivot_wider(id_cols = c(ID_num,session_no,trial),
                  names_from = "method",
                  values_from = "jump_height") |>
      mutate(
        mean = (M2 + M1)/2,
        diff = M2 - M1
      )
    
    lmer_diff <- lmer(diff ~ (1|ID_num/session_no),
                      data = data_wide,
                      REML = TRUE)
    
    totalsd <- sqrt(as.numeric(summary(lmer_diff)$varcor[1])+
                      as.numeric(summary(lmer_diff)$varcor[2])+
                      as.numeric(summary(lmer_diff)$sigma^2)
    )
    
    boot_lmer_diff <- case_bootstrap(lmer_diff, .f = extract_parameters, resample = c(TRUE,FALSE,FALSE), B = 10000)
    
    summary_boot_lmer_diff <- boot_lmer_diff$replicates |>
      mutate(
        totalsd = sqrt(vc1 + vc2 + vc3),
        l_loa = beta - qnorm(1-0.05/2)*totalsd,
        u_loa = beta + qnorm(1-0.05/2)*totalsd
      ) |>
      summarise(
        mean_bias = fixef(lmer_diff)[[1]],
        lower_ci_mean = quantile(beta, 0.025),
        upper_ci_mean = quantile(beta, 0.975),
        lower_ci_l_loa = quantile(l_loa, 0.025),
        upper_ci_l_loa = quantile(l_loa, 0.975),
        lower_ci_u_loa = quantile(u_loa, 0.025),
        upper_ci_u_loa = quantile(u_loa, 0.975),
        l_loa = mean_bias - sqrt(as.numeric(summary(lmer_diff)$varcor[1])+
                                   as.numeric(summary(lmer_diff)$varcor[2])+
                                   as.numeric(summary(lmer_diff)$sigma^2)) * qnorm(1-0.05/2),
        u_loa = mean_bias + sqrt(as.numeric(summary(lmer_diff)$varcor[1])+
                                   as.numeric(summary(lmer_diff)$varcor[2])+
                                   as.numeric(summary(lmer_diff)$sigma^2)) * qnorm(1-0.05/2)
      )
    
    agree_models <- bind_rows(agree_models,
                              tibble(method = i,
                                     summary_boot_lmer_diff))
    
  }
  
  return(agree_models)
  
}


fit_reli_models <- function(data) {
  reli_models <- tibble(method = as.character(),
                        mean_bias = as.numeric(),
                        lower_ci_mean = as.numeric(), 
                        upper_ci_mean = as.numeric(), 
                        lower_ci_l_loa = as.numeric(), 
                        upper_ci_l_loa = as.numeric(),
                        lower_ci_u_loa = as.numeric(),
                        upper_ci_u_loa = as.numeric(),
                        l_loa = as.numeric(),
                        u_loa = as.numeric())
  
  set.seed(1988)
  
  for(i in c("fdi", "jmc", "mj", "opto")) {
    
    data_wide <- data |>
      filter(method == i) |>
      pivot_wider(id_cols = c(ID_num,trial),
                  names_from = "session_no",
                  values_from = "jump_height") |>
      mutate(
        mean = (`2` + `1`)/2,
        diff = `2` - `1`
      )
    
    lmer_diff <- lmer(diff ~ (1|ID_num),
                      data = data_wide,
                      REML = TRUE)
    
    totalsd <- sqrt(as.numeric(summary(lmer_diff)$varcor[1])+
                      as.numeric(summary(lmer_diff)$sigma^2)
    )
    
    boot_lmer_diff <- case_bootstrap(lmer_diff, .f = extract_parameters, resample = c(TRUE,FALSE), B = 10000)
    
    summary_boot_lmer_diff <- boot_lmer_diff$replicates |>
      mutate(
        totalsd = sqrt(vc1 + vc2),
        l_loa = beta - qnorm(1-0.05/2)*totalsd,
        u_loa = beta + qnorm(1-0.05/2)*totalsd
      ) |>
      summarise(
        mean_bias = fixef(lmer_diff)[[1]],
        lower_ci_mean = quantile(beta, 0.025),
        upper_ci_mean = quantile(beta, 0.975),
        lower_ci_l_loa = quantile(l_loa, 0.025),
        upper_ci_l_loa = quantile(l_loa, 0.975),
        lower_ci_u_loa = quantile(u_loa, 0.025),
        upper_ci_u_loa = quantile(u_loa, 0.975),
        l_loa = mean_bias - sqrt(as.numeric(summary(lmer_diff)$varcor[1])+
                                   as.numeric(summary(lmer_diff)$sigma^2)) * qnorm(1-0.05/2),
        u_loa = mean_bias + sqrt(as.numeric(summary(lmer_diff)$varcor[1])+
                                   as.numeric(summary(lmer_diff)$sigma^2)) * qnorm(1-0.05/2)
      )
    
    reli_models <- bind_rows(reli_models,
                             tibble(method = i,
                                    summary_boot_lmer_diff))
    
  }
  
  return(reli_models)
}

##### Funtions for creating results plots

make_agree_plot <- function(data, agree_models) {
  
  plots <- list()
  
  for(i in c("jmc", "mj", "opto")) {
    
    data_wide <- data |>
      filter(method == "fdi" | method == i) |>
      mutate(method = case_when(
        method == "fdi" ~ "M1",
        method == "jmc" ~ "M2",
        method == "mj" ~ "M2",
        method == "opto" ~ "M2"
      )) |>
      pivot_wider(id_cols = c(ID_num,session_no,trial),
                  names_from = "method",
                  values_from = "jump_height") |>
      mutate(
        mean = (M2 + M1)/2,
        diff = M2 - M1
      )
    
    agree_model <- agree_models |>
      filter(method == i)
    
    plot <- data_wide |>
      ggplot(aes(x = mean, y = diff)) +
      
      # Add reference line at zero
      geom_hline(yintercept = 0, linetype = "dashed") +
      
      # Add raw data
      geom_point(alpha = 0.75) +
      
      # Add mean bias
      annotate("rect", alpha = 0.25,
               xmin = -Inf, xmax = Inf,
               ymin = agree_model$lower_ci_mean,
               ymax = agree_model$upper_ci_mean) +
      geom_hline(yintercept = agree_model$mean_bias,
                 size = 1) +
      geom_text(label = glue::glue("Bias = {round(agree_model$mean_bias,2)} cm"),
                x = max(data_wide$mean) + 0.75,
                y = agree_model$upper_ci_mean + 1,
                size = 3,
                # parse = TRUE
      ) +
      geom_text(label = glue::glue("[95%CI: {round(agree_model$lower_ci_mean,2)}, {round(agree_model$upper_ci_mean,2)}]"),
                x = max(data_wide$mean) + 0.75,
                y = agree_model$upper_ci_mean + 0.5,
                size = 3,
                # parse = TRUE
      ) +
      
      # Add lower LoA
      annotate("rect", alpha = 0.25,
               xmin = -Inf, xmax = Inf,
               ymin = agree_model$lower_ci_l_loa[[1]],
               ymax = agree_model$upper_ci_l_loa[[1]]) +
      geom_hline(yintercept = agree_model$l_loa,
                 size = 1, linetype = "dotted")  +
      geom_text(label = glue::glue("Lower LoA = {round(agree_model$l_loa,2)} cm"),
                x = max(data_wide$mean) + 0.75,
                y = agree_model$lower_ci_l_loa - 1,
                size = 3,
                # parse = TRUE
      ) +
      geom_text(label = glue::glue("[95%CI: {round(agree_model$lower_ci_l_loa,2)}, {round(agree_model$upper_ci_l_loa,2)}]"),
                x = max(data_wide$mean) + 0.75,
                y = agree_model$lower_ci_l_loa - 1.5,
                size = 3,
                # parse = TRUE
      ) +
      
      # Add upper LoA
      annotate("rect", alpha = 0.25,
               xmin = -Inf, xmax = Inf,
               ymin = agree_model$lower_ci_u_loa,
               ymax = agree_model$upper_ci_u_loa) +
      geom_hline(yintercept = agree_model$u_loa,
                 size = 1, linetype = "dotted") +
      geom_text(label = glue::glue("Upper LoA = {round(agree_model$u_loa,2)} cm"),
                x = max(data_wide$mean) + 0.75,
                y = agree_model$upper_ci_u_loa + 1,
                size = 3,
                # parse = TRUE
      ) +
      geom_text(label = glue::glue("[95%CI: {round(agree_model$lower_ci_u_loa,2)}, {round(agree_model$upper_ci_u_loa,2)}]"),
                x = max(data_wide$mean) + 0.75,
                y = agree_model$upper_ci_u_loa + 0.5,
                size = 3,
                # parse = TRUE
      ) +
      
      scale_x_continuous(limits = c(min(data_wide$mean), max(data_wide$mean) + 5)) +
      scale_y_continuous(limits = c(min(data_wide$diff)-1, max(data_wide$diff) + 1)) +
      labs(
        x = "Mean of Measurements (cm)",
        y = "Difference of Measurements (cm)"
      ) +
      theme_bw()
    
    
    plots[[i]] <- plot
    
  }
  
  (plots$jmc + labs(title = "Jump Mat")) +
    (plots$mj + labs(title = "MyJump App")) +
    (plots$opto + labs(title = "Optojump")) +
        plot_annotation(title = "Agreement Between Methods",
                    subtitle = "Mixed Effects Model Mean Bias and 95% Limits of Agreement",
                    caption = "Each method compared to gold standard (Force Decks Impulse-Momentum)\n
                    Intervals estimates for mean bias and limits of agreement are 95% quantiles from nonparametric bootstrap resampling participants 10000 times")
  
}

make_reli_plot <- function(data, reli_models) {
  plots <- list()
  
  for(i in c("fdi", "jmc", "mj", "opto")) {
    
    data_wide <- data |>
      filter(method == i)|>
      pivot_wider(id_cols = c(ID_num,trial),
                  names_from = "session_no",
                  values_from = "jump_height") |>
      mutate(
        mean = (`2` + `1`)/2,
        diff = `2` - `1`
      )
    
    reli_model <- reli_models |>
      filter(method == i)
    
    plot <- data_wide |>
      ggplot(aes(x = mean, y = diff)) +
      
      # Add reference line at zero
      geom_hline(yintercept = 0, linetype = "dashed") +
      
      # Add raw data
      geom_point(alpha = 0.75) +
      
      # Add mean bias
      annotate("rect", alpha = 0.25,
               xmin = -Inf, xmax = Inf,
               ymin = reli_model$lower_ci_mean,
               ymax = reli_model$upper_ci_mean) +
      geom_hline(yintercept = reli_model$mean_bias,
                 size = 1) +
      geom_text(label = glue::glue("Bias = {round(reli_model$mean_bias,2)} cm"),
                x = max(data_wide$mean) + 0.75,
                y = reli_model$upper_ci_mean + 1,
                size = 3,
                # parse = TRUE
      ) +
      geom_text(label = glue::glue("[95%CI: {round(reli_model$lower_ci_mean,2)}, {round(reli_model$upper_ci_mean,2)}]"),
                x = max(data_wide$mean) + 0.75,
                y = reli_model$upper_ci_mean + 0.5,
                size = 3,
                # parse = TRUE
      ) +
      
      # Add lower LoA
      annotate("rect", alpha = 0.25,
               xmin = -Inf, xmax = Inf,
               ymin = reli_model$lower_ci_l_loa[[1]],
               ymax = reli_model$upper_ci_l_loa[[1]]) +
      geom_hline(yintercept = reli_model$l_loa,
                 size = 1, linetype = "dotted")  +
      geom_text(label = glue::glue("Lower LoA = {round(reli_model$l_loa,2)} cm"),
                x = max(data_wide$mean) + 0.75,
                y = reli_model$lower_ci_l_loa - 1,
                size = 3,
                # parse = TRUE
      ) +
      geom_text(label = glue::glue("[95%CI: {round(reli_model$lower_ci_l_loa,2)}, {round(reli_model$upper_ci_l_loa,2)}]"),
                x = max(data_wide$mean) + 0.75,
                y = reli_model$lower_ci_l_loa - 1.5,
                size = 3,
                # parse = TRUE
      ) +
      
      # Add upper LoA
      annotate("rect", alpha = 0.25,
               xmin = -Inf, xmax = Inf,
               ymin = reli_model$lower_ci_u_loa,
               ymax = reli_model$upper_ci_u_loa) +
      geom_hline(yintercept = reli_model$u_loa,
                 size = 1, linetype = "dotted") +
      geom_text(label = glue::glue("Upper LoA = {round(reli_model$u_loa,2)} cm"),
                x = max(data_wide$mean) + 0.75,
                y = reli_model$upper_ci_u_loa + 1,
                size = 3,
                # parse = TRUE
      ) +
      geom_text(label = glue::glue("[95%CI: {round(reli_model$lower_ci_u_loa,2)}, {round(reli_model$upper_ci_u_loa,2)}]"),
                x = max(data_wide$mean) + 0.75,
                y = reli_model$upper_ci_u_loa + 0.5,
                size = 3,
                # parse = TRUE
      ) +
      
      scale_x_continuous(limits = c(min(data_wide$mean), max(data_wide$mean) + 5)) +
      scale_y_continuous(limits = c(min(data_wide$diff)-2, max(data_wide$diff) + 2)) +
      labs(
        x = "Mean of Measurements (cm)",
        y = "Difference of Measurements (cm)"
      ) +
      theme_bw()
    
    
    plots[[i]] <- plot
    
  }
  
  (plots$fdi + labs(title = "Force Decks")) +
    (plots$jmc + labs(title = "Jump Mat")) +
    (plots$mj + labs(title = "MyJump App")) +
    (plots$opto + labs(title = "Optojump")) + 
    plot_annotation(title = "Test-Retest Reliability Across Methods",
                    subtitle = "Mixed Effects Model Mean Bias and 95% Limits of Agreement",
                    caption = "Each method compared to gold standard (Force Decks Impulse-Momentum)\n
                    Intervals estimates for mean bias and limits of agreement are 95% quantiles from nonparametric bootstrap resampling participants 10000 times")
  
  
  
} 

make_plot_tiff <- function(plot, path, width, height, device, dpi) {
  
  ggsave(filename = path, plot = plot, width = width, height = height, device = device, dpi = dpi)
  
}
