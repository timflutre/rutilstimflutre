## Contains functions for the "breeding game" Shiny interface
## https://github.com/timflutre/atelier-prediction-genomique

library(shiny)

shinyServer(function(input, output){

  ## tab "plant material"
  readQryPlmat <- reactive({
    if(is.null(input$file.plmat))
      return(NULL)
    read.table(input$file.plmat$datapath, header=TRUE, sep="\t",
               nrows=10)
  })
  output$plmatUploaded <- renderPrint({
    dat <- readQryPlmat()
    if(! is.null(dat))
      ! is.null(readQryPlmat())
  })
  output$plmatSmy <- renderPrint({
    dat <- readQryPlmat()
    if(! is.null(dat))
      summary(dat)
  })
  output$plmatStr <- renderPrint({
    dat <- readQryPlmat()
    if(! is.null(dat))
      str(dat)
  })
  output$qryPlmat <- renderTable({
    readQryPlmat()
  })

  ## tab "phenotypes"
  readQryPheno <- reactive({
    if(is.null(input$file.pheno))
      return(NULL)
    read.table(input$file.pheno$datapath, header=TRUE, sep="\t",
               nrows=10)
  })
  output$phenoUploaded <- renderPrint({
    dat <- readQryPheno()
    if(! is.null(dat))
      ! is.null(readQryPheno())
  })
  output$qryPheno <- renderTable({
    readQryPheno()
  })
  output$dwnlPheno <- downloadHandler(
      filename=function(){
        "phenos.tsv"
      },
      content=function(file){
        write.table(readQryPheno(), file, sep="\t", row.names=FALSE)
      }
  )

  ## tab "genotypes"
  output$qry.geno <- renderTable({
    if(is.null(input$file.geno))
      return(NULL)
    read.table(input$file.geno$datapath, header=TRUE, sep="\t",
               nrows=10)
  })

})
