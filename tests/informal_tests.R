
library(ggplot2)
library(viridis)
library(SimSurvey)
library(plotly)
library(shiny)


sp_pop <- sim_distribution(pop = sim_abundance(),
                           grid = sim_grid(res = c(3.5, 3.5)),
                           space_covar = sim_sp_covar(range = 50, sd = 0.1),
                           ay_covar = sim_ay_covar(sd = 10,
                                                   phi_age = 0.5,
                                                   phi_year = 0.5),
                           depth_par = sim_parabola(alpha = 0, sigma = 50))
# sp_pop <- sim_distribution()
sp_pop$N
d <- sp_pop$sp_N
ages <- sp_pop$ages
years <- sp_pop$years
head(d)

# p <- ggplot(d, aes(x = x, y = y)) +
#   facet_grid(as.numeric(age) ~ as.numeric(year)) +
#   scale_fill_viridis() + theme_void() +
#   theme(panel.spacing = unit(-0.1, "lines"))
# p <- p + geom_tile(aes(fill = prob))
# # p + geom_tile(aes(fill = N))
#
# ggplot(d[sample(seq(nrow(d)), 10000), ]) + geom_point(aes(x = depth, y = N))


x <- sort(unique(d$x))
y <- sort(unique(d$y))
p <- plot_ly(x = x, y = y)

ui <- fluidPage(
  headerPanel("Distribution Explorer"),
  sidebarPanel(
    sliderInput("age", "Age", min = min(ages), max = max(ages),
                value = 1, step = 1),
    sliderInput("year", "Year", min = min(years), max = max(years),
                value = 1, step = 1),
    selectInput("type", "Type", choices = c("surface", "heatmap", "contour"), selected = "surface"),
    checkboxInput("freez", label = "Free z limits?", value = FALSE)
  ),
  mainPanel(
    plotlyOutput("surface", height = "900px", width = "900px")
  )
)

server <- function(input, output) {
  output$surface <- renderPlotly({
    z <- t(xtabs(N ~ x + y, subset = age == input$age & year == input$year, data = d))
    if (input$freez) {
      r <- range(z)
    } else {
      r <- range(d$N)
    }
    add_trace(p, z = z, type = input$type, zauto = input$freez, zmin = min(d$N), zmax = max(d$N)) %>%
      layout(height = input$plotHeight,
             width = input$plotWidth,
             autosize = TRUE,
             scene = list(zaxis = list(range = r)))
  })
}

shinyApp(ui, server)







