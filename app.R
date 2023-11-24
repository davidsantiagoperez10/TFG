## shinydashboard app
library(shiny)
library(shinydashboard)
library(DT)
library(shinycssloaders)
library(shinyWidgets)

# librerías de los métodos
library(haven)
library(tidyr)
library(dplyr)
library(labelled)
library(stringr)
library(FactoMineR)
library(ggplot2)
library(factoextra)
library(rcompanion)
library(ineq)
library(vcd)
library(GoodmanKruskal)
library(corrplot)
library(arules)
library(knitr)
library(hrbrthemes)
library(plotly)
library(forcats)

load("auxiliar.RData")
data <- auxiliar
rm(auxiliar)

### AL INTRODUCIR NUEVOS DATOS: 
# descargar fichero en aplicación (formato .RData), 
# comentar 3 líneas anteriores, y quitarle el comentario a la siguiente:

# data <- readRDS("pr_plays.RData")



# Funciones de análisis

### AC ENTRE 2 VARIABLES
analisis_correspondencias <- function(feature1, feature2){
  
  conting <- table(data[, feature1], data[, feature2])
  conting <- as.matrix(conting)
  
  print(conting)
  cat('\n\n')
  
  # Análisis De Correspondencias (AC)
  res.ca <- CA(conting, graph = TRUE)
  
  print(get_eigenvalue(res.ca))
  
  # Supongamos que en la tabla de contingencias tenemos 26 filas y 8 columnas
  
  # si no hubiera relación, el valor esperado del autovalor para cada 
  # dimensión sería 1/(26-1) = 1/25 = 4% en términos de filas y 1/(8-1) 
  # = 1/7 = 14.286% en términos de columnas.
  
  # Cualquier dimensión que supere el máximo de estos 2 porcentajes 
  # (14.28%) se consideraría importante y útil para sacar conclusiones de 
  # la relación.
  print(fviz_screeplot(res.ca, addlabels = TRUE, ylim = c(0, 100), 
                       main = "% de variabilidad de las dimensiones") +
          geom_hline(yintercept=1/(min(dim(conting))-1)*100, linetype=2, color="red"))
  
  print(fviz_ca_biplot(res.ca, repel = TRUE))
}


### MAPAS DE CALOR: Entre 2 variables

mapas_calor <- function(df = data, col1, col2){
  
  # Creamos la tabla de contingencia de las variables
  aux <- table(df[, col1], df[, col2])
  aux <- as.data.frame(aux)
  
  # Ahora, la transformamos y la representamos con una apariencia 
  # más atractiva a ojos del usuario:
  dfprueba <- aux %>%
    mutate(text = paste0(col1, ": ", Var1, "\n", col2, ": ", Var2, "\n", "Nº Ocurrencias: ", Freq))
  
  p <- ggplot(dfprueba, aes(Var1, Var2, fill= Freq/nrow(df), text=text)) + ggtitle(paste0("Mapa de calor: ", col1, " vs ", col2)) +
    geom_tile() + labs(x = col1, y = col2) +
    geom_text(aes(label = Freq)) + theme_ipsum() +
    scale_fill_gradient('Freq',  low = "gold", high = "red3") 
  
  return(ggplotly(p, tooltip="text"))
}


### MOSAICOS
mosaicos <- function(...){
  cont_table <- table(data[, ...])
  return(mosaicplot(cont_table, main = "Relaciones"))
}


### REGLAS DE ASOCIACIÓN
reglas_asociacion <- function(left_side = ..., right_side = ..., n=10){
  
  # Columnas que analizaremos
  df_tr <- data[, c(left_side, right_side)]
  
  vnum <- which(sapply(data[, c(left_side, right_side)], is.numeric))
  
  # Convertimos las columnas numéricas a factor
  for (i in vnum) {
    df_tr[, names(data[, c(left_side, right_side)])[i]] <- as_factor(df_tr[, names(data[, c(left_side, right_side)])[i]])
  }
  
  # y ahora a un objeto de transacciones
  trans <- as(df_tr, "transactions")
  
  # Determinamos lo que queremos que aparezca en la parte izquierda de la regla
  # (antecedentes) y lo que queremos que salga en la parte derecha (consecuentes)
  aux_lhs <- c()
  for (variable in c(left_side)) {
    aux_lhs <- c(aux_lhs, paste0(variable, "=", unique(df_tr[, variable])))
  }
  print(aux_lhs)
  
  aux_rhs <- c()
  for (variable in c(right_side)) {
    aux_rhs <- c(aux_rhs, paste0(variable, "=", unique(df_tr[, variable])))
  }
  print(aux_rhs)
  
  # Generamos reglas
  rules <- apriori(trans, parameter = list(support = 0.001, confidence = 0.00005), 
                   appearance = list(default = "rhs", lhs = aux_lhs))
  
  # ordenamos y enseñamos
  rules <- sort(rules, by="count", decreasing=TRUE)
  
  # Vemos que las primeras rules tienen la parte izquierda vacía, es decir, 
  # no está considerando ninguna label que nos interesa, por lo que eliminaremos
  # esas reglas, viendo en rules@lhs[i]@data (las reglas son un objeto S4), 
  # que devuelve una matriz que indica qué itemLabel de las reglas se está usando 
  # en la parte izquierda de la transacción, la suma de cada columna, es decir, 
  # el número de labels que usa cada regla en la parte izquierda. Las reglas con 
  # el lhs vacío devolverán 0, por lo que las suprimiremos. En caso de buscar reglas
  # estrictamente entre varias variables sobre una, aumentamos el número de labels
  # mínimo que deban aparecer en la parte izquierda:
  
  
  # índices de las reglas que nos quedaremos
  
  lhs_novacio <- c()
  ########
  for (i in 1:length(rules)) {
    if (sum(rules@lhs[i]@data) == length(left_side)){
      lhs_novacio <- c(lhs_novacio, i)
    }
  }
  
  rules <- rules[lhs_novacio]

  df_rules <- as(rules[1:min(n, length(rules))], "data.frame")
  
  # Mutamos la tabla a una más entendible, ya que sólo aparece una columna
  # del tipo: {var1=x, var2=y => var3=z}, el formato de las reglas en R.
  df_rules <- df_rules %>%
    separate(rules, into = c("lhs", "rhs"), sep = "=>") %>%
    separate(lhs, into = left_side, sep = ",") %>% 
    separate(rhs, into = right_side, sep = ",")
  
  for (col in c(left_side, right_side)) {
    df_rules[, col] <- sub("^.*=([^\\}|\\{]+).*", "\\1", df_rules[, col])
  }
  
  # enseñamos únicamente antecedentes, consecuentes, y conteo de cada regla
  df_rules <- df_rules[, c(left_side, right_side, "count")]
  
  return(datatable(df_rules, 
            options = list(scrollX = "300px", pageLength = 50)))
}


### MCA 
ac_multiple <- function(variables = ...){
  df_red <- data[, variables]
  
  vnum <- which(sapply(df_red, is.numeric))
  
  # Convertimos las columnas numéricas, de existir, a factor
  for (i in vnum) {
    df_red[, names(df_red)[i]] <- as_factor(df_red[, names(df_red)[i]])
  }
  
  res.MCA <- MCA(df_red, graph=FALSE)
  
  # Para que los puntos asociados a las categorías de cada 
  # variable sean distintos por variable
  colores <- c()
  
  for (i in seq_along(variables)) {
    colores <- c(colores, (rep(i, length(levels(df_red[, variables[i]])))))
  }
  
  plot.MCA(res.MCA,invisible= 'ind', col.var=colores,label =c('var'))
}


#-------------------------------


###### APP

ui <- fluidPage(
  list(tags$head(HTML('<link rel="icon", href="https://static.vecteezy.com/system/resources/thumbnails/011/421/099/small/realistic-vector-basketball-isolated-png.png", 
                                   type="image/png" />'))),
  div(style="padding: 10px 0px; width: '30%'",
      titlePanel(
        title="", windowTitle="Análisis P&R"
      )
  ),
  dashboardPage(skin = "black",
    dashboardHeader(title = div(img(src="https://static.vecteezy.com/system/resources/thumbnails/011/421/099/small/realistic-vector-basketball-isolated-png.png", style = "width:5%; height:5%"),
                                "Análisis de jugadas de Pick & Roll"), titleWidth = "40%"),
    dashboardSidebar(
      sidebarMenu(id = "sidebarid",
                  menuItem("Descriptor de Variables", tabName = "var_desc", icon = icon("chart-column"), 
                           actionButton("infodesc", label = "Información", icon = icon("info-circle"), class = "btn-info"),
                           menuSubItem("Descripción individual", tabName = "var_desc_ind"), 
                           menuSubItem("Comparativa con otra variable", tabName = "var_vs_var")),
                  
                  menuItem("Mapas de Calor", tabName = "heatmap", icon = icon("basketball")),
                  conditionalPanel(
                    "input.sidebarid == 'heatmap'",
                    actionButton("info1", label = "Información", icon = icon("info-circle"), class = "btn-info"),
                    
                    selectInput("select_heat", "Variable 1:", choices = names(data)),
                    
                    uiOutput("select_heat_2")
                    ),
                  
                  menuItem("Análisis de Correspondencias", tabName = "ac", icon = icon("basketball")),
                  conditionalPanel(
                    "input.sidebarid == 'ac'",
                    actionButton("info2", label = "Información", icon = icon("info-circle"), class = "btn-info"),
                    
                    selectInput("select_ac", "Variables a estudiar:", choices = names(data), multiple = TRUE, selected = c(names(data)[1], names(data)[2]))
                  ),
                  
                  menuItem("Mosaicos", tabName = "mosaics", icon = icon("basketball")),
                  conditionalPanel(
                    "input.sidebarid == 'mosaics'",
                    actionButton("info4", label = "Información", icon = icon("info-circle"), class = "btn-info"),
                    
                    selectInput("select_mosaic", "Variables a estudiar (hasta 4):", choices = names(data), multiple = TRUE, selected = names(data)[1])
                  ),
                  
                  menuItem("Reglas de asociación", tabName = "arules", icon = icon("basketball"),
                           actionButton("info5", label = "Información", icon = icon("info-circle"), class = "btn-info"),
                           menuSubItem("Reglas", tabName = "arules_rules"),
                           menuSubItem("Contexto", tabName = "arules_context")),
                  conditionalPanel(
                    "input.sidebarid == 'arules_rules'",

                    selectInput("select_arules", "Variables Antecedentes:", choices = names(data), multiple = TRUE, selected = names(data)[1]), 
                    
                    uiOutput("select_arules_2"),
                    
                    numericInput("n", "Nº Reglas a enseñar:", 10, min = 1, max = 100)
                  ), 
                  
                  menuItem("Datos", tabName = "data_full", icon = icon("database"), 
                           menuSubItem("Información", tabName = "data"), 
                           menuSubItem("Añadir más datos", tabName = "data_entry"))
  )),
    dashboardBody(
      tabItems(
        
        tabItem(tabName = "var_desc_ind",
                fluidRow(
                  column(selectInput(inputId = "select_var", label = "Selecciona variable a describir", choices = names(data)), width = 10), 
                  box(
                    h3("Gráfico de frecuencias relativas"),
                    plotOutput("bp_var", height = 300)
                  ), 
                  box(
                    h3("Frecuencias absolutas"),
                    dataTableOutput("levels_var", height = 300)
                  )
                )),
        
        tabItem(tabName = "var_vs_var", 
                fluidRow(
                  column(uiOutput("recordar"), width = 6),
                  column(uiOutput("seg_variable"), width = 6),
                  box(
                    h3("P-Valor del Test Chi-Cuadrado"),
                    verbatimTextOutput("test_chisq"),
                    textOutput("conclusion"),
                    height = 200), 
                  box(
                    h3("V de Cramer"),
                    verbatimTextOutput("test_cramer"),
                    "Esta es la fuerza de la relación entre las variables",
                    height = 200), 
                  box(
                    h3("Gráfico de relación entre categorías"),
                    plotOutput("conting_dispersion"),
                    height = 450), 
                  box(
                    h3("Frecuencias relativas condicionadas"),
                    tableOutput("tabla"),
                    height = 275), 
                  box(
                    h3("Frecuencias absolutas de la segunda variable"),
                    plotOutput("bp_seg"),
                    height = 450)
                )),
  
        tabItem(tabName = "heatmap",
                plotlyOutput("plot_heatmap", height = 750) %>% withSpinner(image = "https://i.gifer.com/origin/e9/e95c7803a86837592c7f8077ea452901.gif")), 
        
        tabItem(tabName = "ac",
                plotOutput("plot_ac", height = 550) %>% withSpinner(image = "https://i.gifer.com/origin/e9/e95c7803a86837592c7f8077ea452901.gif")),
        
        tabItem(tabName = "mosaics",
                plotOutput("plot_mosaic", height = 750) %>% withSpinner(image = "https://i.gifer.com/origin/e9/e95c7803a86837592c7f8077ea452901.gif")),
        
        tabItem(tabName = "arules_rules",
                dataTableOutput("tabla_arules") %>% withSpinner(image = "https://i.gifer.com/origin/e9/e95c7803a86837592c7f8077ea452901.gif")), 
        
        tabItem(tabName = "arules_context",
                fluidRow(
                  column(width=12, h2("Contexto: Frecuencias de las categorías de las variables implicadas")),
                  uiOutput("histograms") %>% withSpinner(image = "https://i.gifer.com/origin/e9/e95c7803a86837592c7f8077ea452901.gif")
                )),

        tabItem(tabName = "data",
                dataTableOutput("tabla_data") %>% withSpinner(image = "https://i.gifer.com/origin/e9/e95c7803a86837592c7f8077ea452901.gif"), 
                fluidRow(
                  box(
                    HTML("<b> Botones de descarga: </b> <br> <br>"),
                    downloadButton("downloadData_csv", "Download Data (.csv)"), 
                    downloadButton("downloadData_rdata", "Download Data (.RData)"),
                  height = 150),
                  box(
                    HTML("<b> Añade más Datos: </b> <br> <br> Diríjase al apartado de <b> Añadir más Datos, </b> y una vez añadidos y visto que el número de filas ha subido, descargue el fichero y cámbielo por el actual, especificado en la línea de <b> load(). </b>"),
                  height = 150)
                )
              ), 
        
        tabItem(tabName = "data_entry",
                fluidRow(
                  box(
                    title = "Jugadas a Añadir",
                    solidHeader = TRUE,
                    width = 12,
                    dataTableOutput("emptyTable")
                  ),
                  box(
                    title = "Acciones",
                    solidHeader = TRUE,
                    width = 12,
                    actionButton("addRow", "Nueva Jugada"),
                    actionButton("submitRows", "Añadir Jugadas"),
                    actionButton("deleteRows", "Borrar")
                  )
                ))
      )
    )
  )
)


server <- function(input, output, session){
  
  #### Botones de información
  
  ## Descriptor de Variables
  observeEvent(input$infodesc, {
    sendSweetAlert(
      session,
      title = "Descriptor de Variables",
      text = HTML(paste("Esta es la sección de", "<b>", "Descriptor de Variables.", "</b>", "<br> <br>", 
                        "<b>", "Estructura:", "</b>", "Presenta dos apartados:  <br> <br>
                        
                        <ul> <li> <b> Descripción individual: </b> Muestra un gráfico de barras de la variable seleccionada, con la frecuencia de cada categoría en la base de datos.
                        A su derecha, tenemos una lista de las categorías en formato de tabla. <br> <br> </li>
                        
                        <li> <b> Comparativa con otra variable: </b> Seleccionada la segunda variable, realizamos las siguientes operaciones: efectuamos el test Chi-Cuadrado de Independencia, 
                        un conocido test estadístico para ver si las variables están asociadas o son independientes; calculamos la V de Cramer para ver cómo de 'fuerte'
                        es la asociación entre las variables, siendo 1 una asociación total, y 0 independencia; mostramos un gráfico de barras apilado en caso de comparar
                        2 variables categóricas, o un diagrama de dispersión en caso de comparar con la variable de resultado; y enseñamos una tabla de proporciones de las
                        categorías de la segunda variable en cada categoría de la primera, o una tabla con el resultado medio de cada categoría. </li> </ul>", "</ul>"
                        
      )),
      type = NULL,
      btn_labels = "Entendido!",
      btn_colors = "blue",
      html = TRUE,
      closeOnClickOutside = TRUE,
      showCloseButton = FALSE,
      width = "60%")
  })
  
  ## Mapas de calor
  observeEvent(input$info1, {
    sendSweetAlert(
      session,
      title = "Mapas de Calor",
      text = HTML(paste("Esta es la sección de", "<b>", "Mapas de Calor.", "</b>", "<br> <br>", 
                        "<ul>", 
                        "<li>", "<b>", "Estructura:", "</b>", "Presenta dos selectores de variables. 
                        Una vez seleccionadas las variables deseadas, se disponen en una cuadrícula donde se enfrentan
                        las categorías de cada una. En cada casilla, se muestra la ocurrencia de ambas categorías a la vez en una jugada.", "</li>", "<br>",

                        "<li>", "<b>", "Explicación:", "</b>", "Este método podría ser útil para ver las tendencias que existen entre
                        categorías en pares de variables. Podríamos extraer conclusiones útiles viendo qué categorías se relacionan en mayor medida con
                        cuáles (la casilla del par de categorías aparecería en rojo), o, caso análogo, qué categorías no se dan a la vez 
                        en una jugada, y, por tanto, no presentan ningún tipo de asociación.
                        ", "</li>", "</ul>", "<br>",
                        
                        "<b>", "¡OJO!: Se recomienda, a la hora de enfrentar 2 variables, ver los gráficos de barras de frecuencias de cada categoría, 
                        ya que es posible que el número de ocurrencias no esté equilibrado, y se saquen conclusiones erróneas debido a categorías que 
                        aparezcan un gran número de veces en las jugadas. Se pueden ver en el", "</b>", "Descriptor de Variables"
      )),
      type = NULL,
      btn_labels = "Entendido!",
      btn_colors = "blue",
      html = TRUE,
      closeOnClickOutside = TRUE,
      showCloseButton = FALSE,
      width = "60%")
  })
  
  ## Análisis de correspondencias
  observeEvent(input$info2, {
    sendSweetAlert(
      session,
      title = "Análisis de Correspondencias",
      text = HTML(paste("Esta es la sección de", "<b>", "Análisis de Correspondencias.", "</b>", "<br> <br>", 
                        "<ul>", 
                        "<li>", "<b>", "Estructura:", "</b>", "Presenta un selector múltiple de variables. 
                        Una vez seleccionadas las variables deseadas, se realiza el análisis, que permitirá explorar
                        y visualizar las relaciones entre las variables gracias al gráfico que se mostrará. Las categorías de
                        variables altamente asociadas aparecerán muy juntas, y alejadas si la relación no es fuerte.", "</li>", "<br>",
                        
                        "<li>", "<b>", "Explicación:", "</b>", "Este método podría verse como una automatización de los mapas de calor, 
                        y una expansión de los mismos al caso multivariante. Analiza las relaciones presentes en la tabla de contingencia (como podríamos ver a un nivel más 'amateur' en los mapas de calor)
                        a través de estadísticos de correspondencia, que calculan la fuerza (cómo de relacionadas están) y dirección (positiva si una categoría provoca que aumente la frecuencia de otra, negativa en caso contrario)
                        de cada relación posible. Posteriormente, se busca una reducción de la dimensionalidad, para poder ver estas relaciones en 2 dimensiones, 
                        recogiendo toda la variabilidad posible.", "</li>", "</ul>", "<br>",
                        
                        "<b>", "¡OJO!: En caso de no estar seleccionada ninguna variable, saldrá un fallo que se soluciona escogiendo variables. En caso de estar seleccionada 1 variable, 
                        se mostrará un gráfico de barras con las frecuencias de sus categorías. Sólo se puede realizar este análisis con 2 o más variables.", "</b>"
      )),
      type = NULL,
      btn_labels = "Entendido!",
      btn_colors = "blue",
      html = TRUE,
      closeOnClickOutside = TRUE,
      showCloseButton = FALSE,
      width = "60%")
  })
  
  ## Mosaicos
  observeEvent(input$info4, {
    sendSweetAlert(
      session,
      title = "Mosaicos",
      text = HTML(paste("Esta es la sección de", "<b>", "Mosaicos.", "</b>", "<br> <br>", 
                        "<ul>", 
                        "<li>", "<b>", "Estructura:", "</b>", "Presenta un selector múltiple de variables. Se podrán escoger hasta 4 variables. En caso de escoger más
                        variables, se cortará el número a las primeras 4. Esto es debido a la pérdida de 'legibilidad' y a la poca capacidad de interpretación que 
                        tienen los gráficos con mayor número de variables. <br> <br>
                        Una vez seleccionadas las variables deseadas, se muestra el mosaico de cuadrículas, que permitirá explorar, para cada par de categorías, la proporción de 
                        valores con respecto al total de cada categoría que presenta una combinación.", "</li>", "<br>",
                        
                        "<li>", "<b>", "Explicación:", "</b>", "Este método permite ver la información que se ve en los mapas de calor, en un concepto de proporciones, y con posibilidad de comparar más de 2 variables.
                        Los pares de categorías cuya cuadrícula sea muy ancha con respecto a las anchuras de su fila, y muy alta con respecto a las alturas de su columna, estarán fuertemente asociadas.", "</li>", "</ul>", "<br>",
                        
                        "<b>", "¡OJO!: En caso de no estar seleccionada ninguna variable, saldrá un fallo que se soluciona escogiendo variables. <br> <br>
                        Se recomienda ver el gráfico de barras de frecuencias de cada categoría para cada variable que se desee estudiar, 
                        ya que es posible que el número de ocurrencias no esté equilibrado, y se saquen conclusiones erróneas debido a categorías que 
                        aparezcan un gran número de veces en las jugadas. Se puede ver cada gráfico en el apartado de", "</b>", "Descriptor de Variables."
      )),
      type = NULL,
      btn_labels = "Entendido!",
      btn_colors = "blue",
      html = TRUE,
      closeOnClickOutside = TRUE,
      showCloseButton = FALSE,
      width = "60%")
  })
  
  ## Reglas de asociación
  observeEvent(input$info5, {
    sendSweetAlert(
      session,
      title = "Reglas de asociación",
      text = HTML(paste("Esta es la sección de", "<b>", "Reglas de asociación.", "</b>", "<br> <br>", 
                        "<ul>", 
                        "<li>", "<b>", "Estructura:", "</b>", "Presenta un selector múltiple y un selector individual de variables, y un selector numérico de número de reglas a enseñar. El selector múltiple se refiere a las variables que deseamos que sean antecedentes, 
                        y el selector individual a la variable que deseamos que sea consecuente. <br> <br>
                        Una vez seleccionadas las variables, se muestra una tabla con los valores de cada variable antecedente y consecuente más frecuente en los datos. Cada fila, que representa una combinación, será una regla de asociación. <br>
                        De esta forma, se podrá ver para qué valores de las variables antecedentes se consigue cada valor de la variable consecuente.", "</li>", "<br>",
                        
                        "<li>", "<b>", "Explicación:", "</b>", "Este método de minería de datos, tal y como está diseñado en la página, permite ver qué combinaciones de variables, llamadas
                        variables antecedentes, dan lugar a cada valor de la variable, generalmente de resultado, llamada consecuente, más veces. Por ejemplo, si queremos ver qué combinación de tipo de defensa, tipo
                        de bloqueo y decisión del jugador con balón se ve más con cada resultado, habría que seleccionarlas como antecedentes, y el resultado como consecuente, y se mostraría una tabla
                        con tantas filas (reglas) como se indique en el selector numérico, con las combinaciones para un resultado determinado con mayor frecuencia.", "</li>", "</ul>", "<br>",
                        
                        "<b>", "¡OJO!: En caso de no estar seleccionada ninguna variable antecedente, saldrá un fallo que se soluciona escogiendo variables. <br> <br>
                        Se recomienda ver el gráfico de barras de frecuencias de cada categoría para cada variable que se desee estudiar, 
                        ya que es posible que el número de ocurrencias no esté equilibrado, y se saquen conclusiones erróneas debido a categorías que 
                        aparezcan un gran número de veces en las jugadas, o, por otra parte, no aparezcan reglas respectivas a categorías con poca frecuencia en los datos 
                        (en ese caso habría que aumenta el número de reglas a enseñar). <br> <br>
                        
                        Para poder tener en cuenta con mayor facilidad este hecho, en el apartado de </b> Contexto <b> se pueden ver los gráficos de barras de cada variable implicada, 
                        y se podrá analizar ese desbalanceo de frecuencias. No obstante, también se puede ver cada gráfico en el apartado de </b>", "Descriptor de Variables."
      )),
      type = NULL,
      btn_labels = "Entendido!",
      btn_colors = "blue",
      html = TRUE,
      closeOnClickOutside = TRUE,
      showCloseButton = FALSE,
      width = "60%")
  })
  
  
  ##########################DESCRIPTOR####################
  
  ## Variable individualmente
  output$bp_var <- renderPlot({
        barplot(table(data[, input$select_var])/nrow(data), xlab=input$select_var, ylab="Frecuencia Relativa", ylim = c(0, 0.9))
  })
  
  
  output$levels_var <- renderDataTable({
    datatable(data.frame(table(data[, input$select_var])), colnames = c("Categoría", "Freq. absoluta"))
  })
  
  
  ## Variable vs otra
  observeEvent(input$select_var, {
    output$recordar <- renderText({
      HTML(paste0("<b>", "Variable Original:", "</b>", "<br>", input$select_var))
    })
  })
  
  observeEvent(input$select_var, {
    output$seg_variable <- renderUI({
      selectInput("seg_variable", label = "Variable con la que comparar:",
                  choices = names(data)[names(data) != input$select_var])
    })
  })
  
  output$test_chisq <- renderPrint({
    chisq.test(table(data[, input$select_var], data[, input$seg_variable]))$p.value
  })
  
  output$conclusion <- renderText({
    if (chisq.test(table(data[, input$select_var], data[, input$seg_variable]))$p.value < 0.05){
      "Podemos rechazar, con un nivel de confianza del 95%, la hipótesis nula que asume independencia entre las variables. 
      Estas variables tienen una asociación significativa."
    } else {
      "Podemos aceptar, con un nivel de confianza del 95%, la hipótesis nula que asume independencia entre las variables"
    }
  })
  
  output$test_cramer <- renderPrint({
    cramerV(table(data[, input$select_var], data[, input$seg_variable]))
  })
  
  output$conting_dispersion <- renderPlot({
    if (is.numeric(data[, input$seg_variable])){
      y.formula <- as.formula(paste(input$seg_variable, input$select_var, sep = "~"))
      boxplot(y.formula, ylab = input$seg_variable, xlab = input$select_var, data = data)
    }
    
    else {
      barplot(100*prop.table(table(data[c(input$seg_variable, input$select_var)]), 2), ylab = "Porcentaje", xlab = input$select_var, legend.text = TRUE)
    }
  })
  
  output$tabla <- renderTable({
    if (is.numeric(data[, input$seg_variable])){
      res <- data[, c(input$select_var, input$seg_variable)] %>% 
        group_by(get(input$select_var)) %>% 
        summarize(mean = mean(get(input$seg_variable)), sd = sd(get(input$seg_variable)))
        
      colnames(res) <- c(input$select_var, "Media", "Desv. Estándar")
      
      res
    
      
    } else {
      as.data.frame.matrix(100*prop.table(table(data[c(input$seg_variable, input$select_var)]), 2))
  }}, rownames = TRUE, colnames = TRUE)
  
  output$bp_seg <- renderPlot({
    barplot(table(data[, input$seg_variable]), xlab=input$seg_variable, ylab="Frecuencia")
  })
  
  
  ##########################1#############################
  ## 1.- HM: elementos de la parte de MAPAS DE CALOR
  observeEvent(input$select_heat, {
    output$select_heat_2 <- renderUI({
        selectInput("select_heat_2", label = "Variable 2:",
                    choices = names(data)[names(data) != input$select_heat])
    })
  })
  
  output$plot_heatmap <- renderPlotly({
    mapas_calor(col1 = input$select_heat, col2 = input$select_heat_2)
  })
  
  ##########################2#############################
  ## 2.- AC: 
  selected <- reactiveValues(values = NULL)
  
  
  observe({
    selected$values <- input$select_ac

    if(length(selected$values) >= 2) {
      
      output$plot_ac <- renderPlot({
        ac_multiple(variables = selected$values)
      })
    }
    
    #else if (length(selected$values) == 2){
     # output$plot_ac <- renderPlot({
      #  analisis_correspondencias(feature1 = selected$values[1], feature2 = selected$values[2])
      #})
    #}
    
    else if (length(selected$values) == 1) {
      output$plot_ac <- renderPlot({
        barplot(table(data[, selected$values]), xlab=names(data)[names(data) == selected$values], ylab="Frecuencia")
      })
    }
    
  })
  
  ##########################3#############################
  ## 3.- MOSAICOS
  selected3 <- reactiveValues(values = NULL)
  
  observe({
    selected3$values <- input$select_mosaic
    
    if(length(selected3$values) > 4) {
      first_4 <- selected3$values[1:4]
      updateSelectInput(session, "select_assoc", selected = first_4)
    }
  })
  
  output$plot_mosaic <- renderPlot({
    mosaicos(input$select_mosaic)
  })
    
  ##########################4#############################
  
  ## 4.- REGLAS de ASOCIACIÓN
  observeEvent(input$select_arules, {
    output$select_arules_2 <- renderUI({
      selectInput("select_arules_2", label = "Variable Consecuente:",
                  choices = names(data)[!(names(data) %in% input$select_arules)],
                  selected = names(data)[names(data) == "resultado"])
    })
  })
  
  # 5.1 Reglas
  output$tabla_arules <- renderDataTable({
    reglas_asociacion(left_side = input$select_arules, right_side = input$select_arules_2, n=input$n)
  })
    
  
  # 5.2 Contexto: histogramas de cada una de las variables
  barplots <- reactiveValues()
  
  observeEvent(input$select_arules, {
    bp_list <- lapply(c(input$select_arules, input$select_arules_2), function(var) {
      ggplot(data, aes_string(x = var)) +
        geom_bar() +
        ggtitle(var)
    })
    
    barplots$plots <- bp_list
  })
    
    output$histograms <- renderUI({
      lapply(seq_along(barplots$plots), function(i) {
        fluidRow(
          column(width = 12, h3(" ")),
          column(width = 12, plotOutput(paste0("plot", i), height = 250, width = "100%"))
        )
      })
    })
    
    observe({
      lapply(seq_along(barplots$plots), function(i) {
        output[[paste0("plot", i)]] <<- renderPlot({
          barplots$plots[[i]]
        })
      })
    })
    

  ##########################5#############################
  
  ## 5.1 Apartado donde se muestran los datos  
  
  datos_orig <- reactiveVal(data)
    
  output$tabla_data <- renderDataTable({
    datatable(datos_orig(), 
              options = list(scrollX = "300px", pageLength = 25))
  })
  
  output$downloadData_csv <- downloadHandler(
    filename = function() {
      paste("pr_plays", ".csv", sep="")
    },
    content = function(file) {
      write.csv(datos_orig(), file)
  })
  
  output$downloadData_rdata <- downloadHandler(
    filename = function() {
      paste("pr_plays", ".RData", sep="")
    },
    content = function(file) {
      saveRDS(datos_orig(), file)
    }, 
    contentType = "application/x-RData")
  
  
  
  ## 5.2 Apartado para añadir datos
  empty_data <- reactiveVal(data.frame(matrix(ncol = ncol(data), nrow = 0, dimnames = list(NULL, colnames(data)))))
  
  # Render the empty datatable
  output$emptyTable <- renderDT({
    datatable(
      empty_data(),
      editable = TRUE,
      options = list(scrollX = "300px", pageLength = 25)
    )
  })
  
  observeEvent(input$addRow, {
    newRow <- reactiveVal(
      lapply(colnames(data), function(col) {
        selectInput(
          inputId = paste0(col, "Input"),
          label = col,
          choices = levels(factor(data[, col]))
        )
      })
    )
    
    showModal(
      modalDialog(
        title = "Nueva Jugada",
        do.call(tagList, newRow()),
        footer = tagList(
          actionButton("saveRow", "Save"),
          modalButton("Cancel")
        )
      )
    )
  })
  
  
  observeEvent(input$saveRow, {
    newValues <- reactiveVal(
      sapply(colnames(data), function(col) {
        input[[paste0(col, "Input")]]
      })
    )
    
    # Append the new row to the existing data
    newData <- rbind(empty_data(), newValues())
    colnames(newData) <- colnames(data)
    
    # Update the reactive dataframe with the new data
    empty_data(newData)
    
    # Close the modal dialog
    removeModal()
  })
  
  
  ## Botón de añadir datos
  observeEvent(input$submitRows, {
    
    ed_tabla <- empty_data()
    
    if (nrow(ed_tabla) > 0){
      newData <- rbind(datos_orig(), ed_tabla)
      
      # Update the data 
      datos_orig(newData)
    
      empty_data(data.frame(matrix(ncol = ncol(data), nrow = 0, dimnames = list(NULL, colnames(data)))))
    }
    
  })
  
  
  observeEvent(input$deleteRows, {
    # Volver a tener el dataframe vacío
    empty_data(data.frame(matrix(ncol = ncol(data), nrow = 0, dimnames = list(NULL, colnames(data)))))
  })
  
  
  
}


shinyApp(ui, server)

