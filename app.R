library(Matrix)
library(DT)
library(shiny)
library(shinybusy)
library(Seurat)
library(dplyr)
library(shinyWidgets)
library(ggplot2)
library(scater)

### set a size limit for uploading files
options(shiny.maxRequestSize = 3000*1024^2)
### set a seed
set.seed(42)

### UI
{
  ui <- fixedPage(
    
    
    fixedRow(
      titlePanel(paste("visualisation of snRNA-seq data set"))
    ),
    hr(),
    hr(),
    
    ### load the data
    
    fixedRow(
      titlePanel("upload the data"),
      HTML(paste("All data should be in .rda format.", 
        "1) Seurat object should have PCA and UMAP slots, data should be clustered." ,
        "2) DE results should contain a list of tables - one per clusters, generated with Seurat::FindMarkers() function." , 
        "3) Cluster markers file should contain one table generated with Seurat::FindAllMarkers", sep="<br/>")) ),
    h1(),
    
    fixedRow(
      column(8,fileInput("uploadSeurat", NULL, width = "70%", buttonLabel = "Upload Seurat object", multiple = F)),
      column(8,fileInput("uploadDEresults", NULL, width = "70%", buttonLabel = "Upload results of DE analysis", multiple = F)),
      column(8,fileInput("uploadClustMarks", NULL, width = "70%", buttonLabel = "Upload cluster markers", multiple = F))
      ),
    fixedRow( column(4,actionButton( inputId="do","Click to plot")),
              add_busy_spinner(spin = "fading-circle") #progress indicator
    ),
    hr(),
    hr(),
    
    ### UMAPS
    fixedRow(
      titlePanel("Dataset and Cluster Metadata Inspection"),
      p(paste(
        "A look at the whole dataset.",
        "The dimensional reduction is UMAP (cells",
        "projected into 2D space), where proximity indicates transcriptional similarity.\n",
        "With the option on top you can add a metadata overlay",
        "to the cell projection. "
      )),
   
      h1()
    ),
    
    hr(),
  
    fixedRow(
      column(4,uiOutput("umapMDcol")),
      # column(2,uiOutput("integratedVIZ")),
      column(2,uiOutput("umapMDlog"))
    ),
    
    fixedRow(
      column(10,plotOutput("umap",height="700px"))
    ),
    
    fixedRow(
      column(10,align="right",downloadButton("umapSave",paste("Save image")))
      ) ,
    hr(),
    hr(),
    
    ### scatter plots
    fixedRow(
      p(paste("plot continuous meta data features"
      )),
      h1()),
    
    fixedRow(
      column(2,uiOutput("mdScatterX")),
      column(2,uiOutput("mdScatterY")), 
      column(2,uiOutput("scatterLog"))
    ),
    fixedRow(
      column(12,plotOutput("mdScatter",height="500px"))
    ),
    fixedRow(
      column(12,align="right",downloadButton("mdScatterSave",paste("Save as pdf"))),
    ) ,
    
    hr(),
    hr(),
    
    ### violin plots of gene expression
    fixedRow(
      titlePanel("Violin plots of Gene Expression"),h1()
    ) ,
    
    hr(),
    
    fixedRow(
      column(3,searchInput(
        inputId = "GeneList", 
        label = "Enter list of genes for violin plot(comma-separated):", 
        value = c("Wt1,Epha6"), 
        btnSearch = icon("search"), 
        btnReset = icon("remove"), 
        width = "100%"
      ))
      ,
      column(2,uiOutput("ListType")),
      column(2,uiOutput("showcells"))
      ),
    
    fixedRow(
      column(4,uiOutput("AllorOne1")),
      column(2,uiOutput("WhichCluster1"))
    ),
    
    fixedRow(
      column(12,plotOutput("vplot",height="600px"))
    ),
    
    fixedRow(
      column(12,align="right",downloadButton("vSave",paste("Save image"))))
    ,
    
    hr(),
    hr(),
    
    ### gene expression lvls in individual cells
    fixedRow(
      titlePanel("Expression of Genes of Interest across individual cells"),h1()
    ) ,
    hr(),
    
    fixedRow(
      column(3,searchInput(
        inputId = "GeneName", 
        label = "Enter gene name :", 
        value = "Wt1", 
        btnSearch = icon("search"), 
        btnReset = icon("remove"), 
        width = "100%"
      ))
    ) ,
    
    fixedRow(
      column(4,uiOutput("AllorOne2")),
      column(2,uiOutput("WhichCluster2"))
    ),
    
    fixedRow(
      column(6,plotOutput("goiplot1",height="580px")),
      column(6,plotOutput("goiplot2",height="580px"))
    ),
    
    
    fixedRow(
      column(12,align="right",downloadButton("geneSave",paste("Save image")))
      ) ,
    hr(),
    hr(),
    
    ### cluster markers UI
    fixedRow(
      titlePanel("Cluster Markers")
    ) ,
    hr(),
    
    fixedRow(
      p("Expression heatmap of cluster markers, top 10 markers (max) per cluster are shown")
    ), 
    
    fixedRow(
      column(12,plotOutput("clustPlot",height="1400px"))
    ) ,
    fixedRow(
      column(12,align="right",downloadButton("hmapSave",paste("Save image")))
    ) ,
    
    hr(),
   basicPage(
      DT::dataTableOutput("mytable")
    ) ,
   h1(),
    hr(),
    hr(),

    
   
    
    ### UI for DE  analysis
    fixedRow(
      titlePanel("Differential Expression between various groups of samples")
      # p("The last column shows what groups were tested against each other")
    ) ,
    hr(),
    
    basicPage(
      DT::dataTableOutput("DEtable")
    ) ,
    hr(),
    hr()
   
  
  )
}


####=========================SERVER=====================####
server <- function(input, output) {
  
  ## load the Seurat data
  sce.combined <- reactive({
    req( input$uploadSeurat )
    readRDS( file=input$uploadSeurat$datapath ) 
  })

  
  ## load the DE results
  sce_DEresults <- reactive({
    req( input$uploadDEresults )
    readRDS( file=input$uploadDEresults$datapath ) 
  })

  
  ## load the cluster markers results
  sce_ClustMarks <- reactive({
    req( input$uploadClustMarks )
    readRDS( file=input$uploadClustMarks$datapath ) 
  })
  
  
 # start rendering only after clicking the button
  observeEvent( input$do, {
    
    show_spinner()
    
  ### UMAPs
  {
      output$umapMDcol <- renderUI({
        sce.combined <- sce.combined()
        selectInput("umapMDcol",label="Metadata:",width="60%",choices=colnames( sce.combined@meta.data),
                    selected=colnames( sce.combined@meta.data)[8])
      })
      
      output$umapMDlog <- renderUI({
        selectInput("umapMDlog",label="Type of reduction:",width="60%",choices=list("UMAP"="umap","PCA"="pca"),
                    selected="UMAP")
      })
      
      output$umap <- renderPlot({
        req( input$umapMDcol )
        req( input$umapMDlog )
        
        sce.combined <- sce.combined()
        
        if ( is.character( sce.combined@meta.data[ ,input$umapMDcol]) | is.factor( sce.combined@meta.data[ ,input$umapMDcol]) ) {
          print( DimPlot( sce.combined, pt.size=1, reduction=input$umapMDlog, group.by= input$umapMDcol , label=T ,
                          shuffle = T ))
        } else {
          print(FeaturePlot( sce.combined, pt.size=1, label=F, reduction=input$umapMDlog,
                             features= input$umapMDcol,  cols=c("lightgrey", "darkblue")))}
        
      })
      
      output$umapSave <- downloadHandler(filename=paste0("umap",input$umapMDcol,".pdf"),
                                         content=function(file){
                                           req(input$umapMDcol)
                                           
                                           sce.combined <- sce.combined()
                                          
                                           grDevices::cairo_pdf(file,height=10,width=10,fallback_resolution=1200)


                                           if ( is.character( sce.combined@meta.data[ ,input$umapMDcol]) | is.factor( sce.combined@meta.data[ ,input$umapMDcol]) ) {
                                             print( DimPlot( sce.combined, pt.size=1, label=F, reduction=input$umapMDlog,
                                                             group.by= input$umapMDcol  , shuffle = T ))
                                           } else {
                                             print(FeaturePlot( sce.combined, pt.size=1, label=F, reduction=input$umapMDlog,
                                                                features= input$umapMDcol , cols=c("lightgrey", "darkblue") ))}

                                           grDevices::dev.off()
                                         })
    }
    

  ### scatter plots of continuous meta-features
  {
    output$mdScatterX <- renderUI({
      # update the dataset
      sce.combined <- sce.combined()
      
      selectInput(
        "mdScatterX","X axis:",choices=colnames( sce.combined@meta.data)[c(2,3,7,9,10)],
        selected=colnames( sce.combined@meta.data)[2]
      )
    })

    output$mdScatterY <- renderUI({
      # update the dataset
      sce.combined <- sce.combined()
      
      selectInput(
        "mdScatterY","Y axis:",choices=colnames( sce.combined@meta.data)[c(2,3,7,9,10)],
        selected=colnames( sce.combined@meta.data)[3]
      )
    })

    output$scatterLog <- renderUI({
      temp_choices <- list("Yes"="Yes","No"="No")


      radioButtons("scatterLog",inline=F,label="Group by original identity?",
                   choices= temp_choices, selected="No")

    })

    output$mdScatter <- renderPlot({
      req(input$scatterLog)
      
      # update the dataset
      sce.combined <- sce.combined()
      
      if(input$scatterLog=="No"){
        pp1 <- (FeatureScatter( sce.combined, feature1=input$mdScatterX,feature2=input$mdScatterY, group.by=NULL))
      }
      if(input$scatterLog=="Yes"){
        pp1 <- (FeatureScatter( sce.combined, feature1=input$mdScatterX,feature2=input$mdScatterY, group.by="orig.ident"))
      }
      pp2 <- ggplot( data= sce.combined@meta.data[order( sce.combined@meta.data[,input$mdScatterX]), ] ,
                     aes(y=  sce.combined@meta.data[order( sce.combined@meta.data[,input$mdScatterX]),input$mdScatterX], x=(1:nrow( sce.combined@meta.data))) ) +
        geom_point() + theme_bw() + xlab("index")+ ylab( input$mdScatterX )
      pp3 <- ggplot( data=  sce.combined@meta.data[order( sce.combined@meta.data[,input$mdScatterY]), ] ,
                     aes(y=  sce.combined@meta.data[order( sce.combined@meta.data[,input$mdScatterY]),input$mdScatterY],
                         x=(1:nrow( sce.combined@meta.data))) ) +
        geom_point() + theme_bw() + xlab("index")+ ylab( input$mdScatterY)
      print( cowplot::plot_grid( pp2, pp3, pp1 , nrow = 1))
    })

    output$mdScatterSave <- downloadHandler( filename=paste0("mdScatter.pdf"),
                                             content=function(file){
                                               
                                               req(input$scatterLog)
                                               # update the dataset
                                               sce.combined <- sce.combined()
                                               
                                               if(input$scatterLog=="No"){
                                                 pp1 <- (FeatureScatter( sce.combined, feature1=input$mdScatterX,feature2=input$mdScatterY, group.by=NULL))
                                               }
                                               if(input$scatterLog=="Yes"){
                                                 pp1 <- (FeatureScatter( sce.combined, feature1=input$mdScatterX,feature2=input$mdScatterY, group.by="orig.ident"))
                                               }
                                               pp2 <- ggplot( data= sce.combined@meta.data[order( sce.combined@meta.data[,input$mdScatterX]), ] ,
                                                              aes(y=  sce.combined@meta.data[order( sce.combined@meta.data[,input$mdScatterX]),input$mdScatterX], x=(1:nrow( sce.combined@meta.data))) ) +
                                                 geom_point() + theme_bw() + xlab("index")+ ylab( input$mdScatterX )
                                               pp3 <- ggplot( data=  sce.combined@meta.data[order( sce.combined@meta.data[,input$mdScatterY]), ] ,
                                                              aes(y=  sce.combined@meta.data[order( sce.combined@meta.data[,input$mdScatterY]),input$mdScatterY],
                                                                  x=(1:nrow( sce.combined@meta.data))) ) +
                                                 geom_point() + theme_bw() + xlab("index")+ ylab( input$mdScatterY)

                                               grDevices::cairo_pdf(file,height=6,width=12,fallback_resolution=1200)
                                               print( cowplot::plot_grid( pp2, pp3, pp1 , nrow = 1))
                                               grDevices::dev.off()
                                             }  )
  }

  ### violin plots of genes of interest
  {
    output$AllorOne1 <- renderUI({
      selectInput("AllorOne1",label="Markers of all clusters or one?",width="100%",choices=list("All"="all2","One"="one2"),
                  selected="all2")
    })

    output$WhichCluster1 <- renderUI({
      req(input$AllorOne1)
      
      # update the dataset
      sce.combined <- sce.combined()
      
      if(input$AllorOne1=="one2"){
        selectInput("WhichCluster1",label="Which Cluster?", width="100%",choices=levels(  sce.combined$seurat_clusters),
                    selected="0", multiple = T)
      }else{return(NULL)}
    })

    output$ListType <- renderUI({
      temp_choices <- list("Clusters"="clusterlist", "Mouse genotype"="gtypelist")

      radioButtons("ListType",inline=F,label="What kind of groups do you want to use?",
                   choices=temp_choices,selected="clusterlist")

    })

    output$showcells <- renderUI({
      temp_choices <- list("Yes"="Yes","No"="No")

      radioButtons("showcells",inline=F,label="show individual cells",
                   choices= temp_choices, selected="No")
    })

    output$vplot <- renderPlot({
      req(input$ListType)
      
      # split gene names by comma
      translate <- unlist( strsplit( input$GeneList , split = ",") )

      # update the dataset
      sce.combined <- sce.combined()
      
      ## all or selected clusters
      if ( input$AllorOne1=="all2") datTOplot <-  sce.combined else {
        datTOplot <- subset(  sce.combined , ident= input$WhichCluster1  )
      }
      ## show all cells or not
      if ( input$showcells=="No") pointsize	<- 0 else pointsize	<- 0.1

      ## type of groupping
      if(input$ListType=="clusterlist"){
        print(VlnPlot(object = datTOplot , assay =  "RNA", features = translate, pt.size=pointsize, group.by="seurat_clusters") )
      }
            if(input$ListType=="origidentlist"){
        print(VlnPlot(object = datTOplot, assay =  "RNA", features = translate, pt.size=pointsize,group.by="orig.ident"))
      }
      if(input$ListType=="gtypelist"){
        print(VlnPlot(object = datTOplot, assay =  "RNA", features = translate, pt.size=pointsize,group.by="gtype"))
      }
    })

    output$vSave <- downloadHandler(filename=paste0("vplot.pdf"),
                                    content=function(file){
                                      require(input$GeneList)
                                      grDevices::cairo_pdf(file,height=5,width=10,fallback_resolution=1200)
                                      req(input$ListType)
        
                                      translate <-unlist( strsplit( input$GeneList , split = ",") )

                                      # update the dataset
                                      sce.combined <- sce.combined()
                                      
                                      if ( input$AllorOne1=="all2") datTOplot <-  sce.combined else {
                                        datTOplot <- subset(  sce.combined , ident= input$WhichCluster1 )
                                      }

                                      ## show all cells or not
                                      if ( input$showcells=="No") pointsize	<- 0 else pointsize	<- 0.1

                                      ## type of groupping
                                      if(input$ListType=="clusterlist"){
                                        print(VlnPlot(object = datTOplot , assay =  "RNA", features = translate, pt.size=pointsize, group.by="seurat_clusters") )
                                      }
                                      if(input$ListType=="origidentlist"){
                                        print(VlnPlot(object = datTOplot, assay =  "RNA", features = translate, pt.size=pointsize,group.by="orig.ident"))
                                      }
                                      if(input$ListType=="gtypelist"){
                                        print(VlnPlot(object = datTOplot, assay =  "RNA", features = translate, pt.size=pointsize,group.by="gtype"))
                                      }
                                      grDevices::dev.off()
                                    } )
  }

  ###  expression of genes of interest in individual cells
  {

    # UMAP with cell-type labeles stored in cell_type metadata column
    output$goiplot1 <- renderPlot({
      
      # update the dataset
      sce.combined <- sce.combined()
      
      # make a plot
      plot3<- DimPlot(  sce.combined, group.by="seurat_clusters" )
      
      # add celltype labels
      plot3$data$Label <-  sce.combined$cell_type
      print(LabelClusters(plot3, id="Label"))
      
    })

    # choose if all cells or one cluster should be vizualised
    output$AllorOne2 <- renderUI({
      selectInput("AllorOne2",label="Display all clusters or one?",width="100%",choices=list("All"="all3","One"="one3"),
                  selected="all3")
    })
    # if one cluster - which one?
    output$WhichCluster2 <- renderUI({
      req(input$AllorOne2)
      
      # update the dataset
      sce.combined <- sce.combined()
      
      if(input$AllorOne2=="one3"){
        selectInput("WhichCluster2",label="Which Cluster?", width="100%",choices=levels(  sce.combined$seurat_clusters ),
                    selected="0", multiple = T)
      }else{return(NULL)}
    })
    # UMAP showing expression of individual genes
    output$goiplot2 <- renderPlot({

      # update the dataset
      sce.combined <- sce.combined()
      # set default to RNA
      DefaultAssay(sce.combined) <- "RNA"

      if ( input$AllorOne2 =="all3") {
        datTOplot <-  sce.combined
        ppsize=1} else {
        datTOplot <- subset(  sce.combined , ident= input$WhichCluster2 )
        ppsize=3
      }

      print( scater::plotUMAP( as.SingleCellExperiment( datTOplot), colour_by= input$GeneName, point_size = ppsize   ) )
    })

    # save both UMAPs
    output$geneSave <- downloadHandler(filename=paste0(input$GeneName,".pdf"),
                                       content=function(file){
                                         require(input$GeneName)
                                         
                                         # update the dataset
                                         sce.combined <- sce.combined()
                                         
                                         DefaultAssay(sce.combined) <- "RNA"

                                         # all or selected clusters
                                         if ( input$AllorOne2 =="all3") {
                                           datTOplot <-  sce.combined
                                           ppsize=1 } else {
                                             ppsize=4
                                           datTOplot <- subset(  sce.combined , ident= input$WhichCluster2 )
                                         }

                                         grDevices::cairo_pdf(file,height=6,width=12,fallback_resolution=1200)

                                         plot3<- DimPlot(  sce.combined, group.by="seurat_clusters")
                                         plot3$data$Label <-  sce.combined$seurat_clusters
                                         pp1 <- ( LabelClusters(plot3, id="Label"))
                                         pp2 <- ( scater::plotUMAP( as.SingleCellExperiment( datTOplot),
                                                                    colour_by= input$GeneName, point_size = ppsize   ) )
                                         print( cowplot::plot_grid( pp1,pp2) )
                                         grDevices::dev.off()
                                       }  )
  }

  ### plot and table with cluster markers
  {
    # render heatmap of cluster markers
    output$clustPlot <- renderPlot({
      # update the data
      sce_ClustMarks <- sce_ClustMarks()
      sce.combined <- sce.combined()
      
      # select top 10 markers for each cluster
      top10 <- sce_ClustMarks %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
      
      # subset the data
      sce.combined_subs <- subset( sce.combined ,  DoubletScore< 1 , downsample = 200  )
      sce.combined_subs <- ScaleData( sce.combined_subs, features = rownames( sce.combined_subs) )
      
      # make a plot
      print( DoHeatmap(  sce.combined_subs, features = top10$gene) + NoLegend() ) 
      })
    # save button for the heatmap of cluster markers
    output$hmapSave <- downloadHandler( filename=paste0("clusterMarks_heatmap.pdf"),
                                       content=function(file){
                                         
                                         # update the data
                                         sce_ClustMarks <- sce_ClustMarks()
                                         sce.combined <- sce.combined()
                                         # select top 10 markers for each cluster
                                         top10 <- sce_ClustMarks %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
                                         # subset the data
                                         sce.combined_subs <- subset( sce.combined ,  DoubletScore< 1 , downsample = 200  )
                                         sce.combined_subs <- ScaleData( sce.combined_subs, features = rownames( sce.combined_subs) )
                                         
                                         
                                       grDevices::cairo_pdf(file,height=18,width=15,fallback_resolution=1200)
                                         print( DoHeatmap(  sce.combined_subs, features = top10$gene) + NoLegend() )
                                         grDevices::dev.off()
                                       }  )
    # render table of cluster markers
    output$mytable <- DT::renderDataTable({

      # update the data
      sce_ClustMarks <- sce_ClustMarks()
      # visualise
      datatable( sce_ClustMarks , rownames = F )   } , 
        extensions= 'Buttons', options = list( scrollY = 500 , dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),paging=F))

  }

  ### DE table
  {
    # render table of DE for all clusters
    output$DEtable <- DT::renderDataTable({

      # update the data
      sce_DEresults <- sce_DEresults()
      sce.combined <- sce.combined()
      
      # make a tab with all clusters combined
      allTab <- Reduce( rbind, lapply( seq(sce_DEresults), function(ii) {
          sce_DEresults[[ii]]$cluster <- levels(  sce.combined$seurat_clusters )[ii]
          return(sce_DEresults[[ii]])
        }) )

      names(sce_DEresults) <- levels(  sce.combined$seurat_clusters )

      
      datatable( allTab, rownames = F )  } ,  extensions= 'Buttons', options = list( scrollY = 500 , dom = 'Bfrtip',
                                                buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),paging=F))


  }
  

} )

}

shinyApp(ui = ui, server = server)



