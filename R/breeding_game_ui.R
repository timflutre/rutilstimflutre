## Contains functions for the "breeding game" Shiny interface
## https://github.com/timflutre/atelier-prediction-genomique

## test locally:
## library(shiny)
## setwd("<...>") # path to dir containing a dir named "breeding-game"
## runApp("breeding-game") # make symlinks ui.R and server.R

## or when distributed in a package:
## shiny::runApp(system.file('appdir', package='packagename'))

library(shiny)

shinyUI(
    navbarPage(
        title="Joue au sélectionneur !",

        tabPanel(title="Introduction",
                 sidebarLayout(
                     sidebarPanel(
                         p("Viens jouer à sélectionner de nouvelles variétés de plantes !")
                     ),
                     mainPanel(
                         h1("Informations biologiques"),
                         h2("Apimeta simulans, une espèce pleine d’avenir"),
                         p("Découverte récemment aux confins de la vallée supérieure de l’Aghromonpe, Apimeta simulans appartient au genre des Statisticeae. Elle produit des fleurs qui contiennent un composé alcaloïde, la sepmetine, consommée par les étudiants pour éviter les maux de tête lors des efforts intellectuels trop intenses. Le marché est donc très important et se développe rapidement. Apimeta peut être sensible à quelques maladies, dont la redoutable rouille fluorescente, Putrida psychedelica. Les producteurs sont payés à la quantité produite (les rendements sont de l’ordre de 40 kg de fleurs par ha), mais les transformateurs ont réussi à exiger que la teneur moyenne en sepmetine des lots commerciaux soit au dessus de 14 pour mille."),
                         h2("Biologie de l’espèce"),
                         p("Apimeta simulans est hermaphrodite et autogame. Elle est produite sous forme de lignées pures, et se croise facilement. On peut produire en serre jusqu’à deux générations par an pour accélerer la fixation et le retour à l’homozygotie. Il est également possible de produire des haploides doublés. Le taux de multiplication de l’espèce est très élevé, chaque plante étant capable de produire plus de 1000 graines."),
                         h2("...")
                     )
                 )),

        tabPanel(
            title="Croiser",
            sidebarLayout(
                sidebarPanel(
                    fileInput(inputId="file.plmat",
                              label=h3("Choisis un fichier:"),
                              multiple=FALSE,
                              accept=c(".txt", ".tsv"))
                ),
                mainPanel(
                    tabsetPanel(type="tabs",
                                tabPanel("vérif",
                                         verbatimTextOutput("plmatUploaded")),
                                tabPanel("résumé",
                                         verbatimTextOutput("plmatSmy"),
                                         verbatimTextOutput("plmatStr")),
                                tabPanel("données",
                                         tableOutput(outputId="qryPlmat")),
                                tabPanel("sortie", p("RAS"))
                                )
                )
            )),

        tabPanel(
            title="Phénotyper",
            sidebarLayout(
                sidebarPanel(
                    fileInput(inputId="file.pheno",
                              label=h3("Choisis un fichier:"),
                              multiple=FALSE,
                              accept=c(".txt", ".tsv"))
                ),
                mainPanel(
                    tabsetPanel(
                        tabPanel("vérif",
                                 verbatimTextOutput("phenoUploaded")),
                        tabPanel("résumé",
                                 p("TODO")),
                        tabPanel("données", tableOutput(outputId="qryPheno")),
                        tabPanel("sortie", downloadButton("dwnlPheno",
                                                          "Télécharger")))
                )
            )),

        tabPanel(
            title="Génotyper",
            sidebarLayout(
                sidebarPanel(
                    fileInput(inputId="file.geno",
                              label=h3("Choisis un fichier:"),
                              multiple=FALSE,
                              accept=c(".txt", ".tsv"))
                ),
                mainPanel(
                    p("TODO")
                )
            )),

        tabPanel(title="A propos",
                 sidebarLayout(
                     sidebarPanel(
                         p("Auteurs: Timothée Flutre, Jacques David"),
                         p("Copyright 2014-2017 INRA, Montpellier SupAgro")
                     ),
                     mainPanel(
                         p("L’accès au génome permet de produire de nombreuses données de génotypage sur un nombre important d’individus qui peuvent être phénotypés à des caractères d'intérêt agronomique. L’association entre marqueurs et caractères fait souvent l’impasse sur les gènes ayant de faibles effets, ceux-ci pouvant pourtant déterminer une part importante de la variation d'origine génétique, et donc du  progrès génétique envisageable. La prédiction génomique des valeurs des plantes, puis l’utilisation de ces prédictions dans un programme de sélection, est donc un des enjeux majeurs."),
                         br(),
                         p("Il n’est pas facile d’enseigner ces nouveaux aspects car, d’une part, d'importants pré-requis sont nécessaires, et d’autre part, il n’existe pas de méthodologie définitive valable pour tous modèles biologiques, avec des enjeux variables en terme de programmation de progrès génétique et de sortie variétale. C’est un domaine très actif de recherche avec de nombreuses innovations en statistique et génétique."),
                         p(" Notre ambition est de confronter les étudiants à certains aspects théoriques et pratiques dans une dynamique d’apprentissage actif basée sur la prise en main pratique d’outils et sur l’utilisation d’un \"jeu sérieux\"."),
                         br(),
                         p("Les documents pédagogiques sont disponibles sur ",
                           a("GitHub.",
                             href="https://github.com/timflutre/atelier-prediction-genomique"))
                     )
                 ))
    )
)
