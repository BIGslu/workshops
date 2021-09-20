library(kimma)
#Linear models
model_result <- kmFit(dat = example.voom, kin = example.kin, patientID = "donorID",
                      run.lmekin = TRUE, run.lme = TRUE,
                      model = "~ virus + (1|donorID)", processors = 6)

#plot fxn
sigma_plot <- function(model_result, x, y){
  library(tidyverse)

  #Extract results
  dat_x <- model_result[[x]] %>%
    distinct(model, gene, sigma)
  dat_y <- model_result[[y]] %>%
    distinct(model, gene, sigma)
  #Merge and format
  dat <- bind_rows(dat_x,dat_y) %>%
    pivot_wider(names_from = model, values_from = sigma) %>%
    #add best fit variable
    mutate(`Best fit` = ifelse(get(x)<get(y), x,
                               ifelse(get(y)<get(x), y, "none")))

  #plot
  plot <- ggplot(dat, aes(x=get(x), y=get(x), color=`Best fit`)) +
    geom_point(alpha=0.3) +
    labs(x=x, y=y) +
    geom_abline() +
    theme_classic()


  message("Total genes best fit by")
  table(dat$`Best fit`)
  return(plot)
}

#test fxn
sigma_plot(model_result, x="lme", y="lmekin")
