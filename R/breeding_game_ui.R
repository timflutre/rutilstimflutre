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
        title="Joue au s\u00e9lectionneur !",

        tabPanel(title="Introduction",
                 sidebarLayout(
                     sidebarPanel(
                         p("Viens jouer \u00e0 s\u00e9lectionner de nouvelles vari\u00e9t\u00e9s de plantes !")
                     ),
                     mainPanel(
                         h1("Informations biologiques"),
                         h2("Apimeta simulans, une esp\u00e8ce pleine d'avenir"),
                         p("D\u00e9couverte r\u00e9cemment aux confins de la vall\u00e9e sup\u00e9rieure de l'Aghromonpe, Apimeta simulans appartient au genre des Statisticeae. Elle produit des fleurs qui contiennent un compos\u00e9 alcalo\u00efde, la sepmetine, consomm\u00e9e par les \u00e9tudiants pour \u00e9viter les maux de t\u00eate lors des efforts intellectuels trop intenses. Le march\u00e9 est donc tr\u00e8s important et se d\u00e9veloppe rapidement. Apimeta peut \u00eatre sensible \u00e0 quelques maladies, dont la redoutable rouille fluorescente, Putrida psychedelica. Les producteurs sont pay\u00e9s \u00e0 la quantit\u00e9 produite (les rendements sont de l'ordre de 40 kg de fleurs par ha), mais les transformateurs ont r\u00e9ussi \u00e0 exiger que la teneur moyenne en sepmetine des lots commerciaux soit au dessus de 14 pour mille."),
                         h2("Biologie de l'esp\u00e8ce"),
                         p("Apimeta simulans est hermaphrodite et autogame. Elle est produite sous forme de lign\u00e9es pures, et se croise facilement. On peut produire en serre jusqu'\u00e0 deux g\u00e9n\u00e9rations par an pour acc\u00e9lerer la fixation et le retour \u00e0 l'homozygotie. Il est \u00e9galement possible de produire des haploides doubl\u00e9s. Le taux de multiplication de l'esp\u00e8ce est tr\u00e8s \u00e9lev\u00e9, chaque plante \u00e9tant capable de produire plus de 1000 graines."),
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
                                tabPanel("v\u00e9rif",
                                         verbatimTextOutput("plmatUploaded")),
                                tabPanel("r\u00e9sum\u00e9",
                                         verbatimTextOutput("plmatSmy"),
                                         verbatimTextOutput("plmatStr")),
                                tabPanel("donn\u00e9es",
                                         tableOutput(outputId="qryPlmat")),
                                tabPanel("sortie", p("RAS"))
                                )
                )
            )),

        tabPanel(
            title="Ph\u00e9notyper",
            sidebarLayout(
                sidebarPanel(
                    fileInput(inputId="file.pheno",
                              label=h3("Choisis un fichier:"),
                              multiple=FALSE,
                              accept=c(".txt", ".tsv"))
                ),
                mainPanel(
                    tabsetPanel(
                        tabPanel("v\u00e9rif",
                                 verbatimTextOutput("phenoUploaded")),
                        tabPanel("r\u00e9sum\u00e9",
                                 p("TODO")),
                        tabPanel("donn\u00e9es", tableOutput(outputId="qryPheno")),
                        tabPanel("sortie", downloadButton("dwnlPheno",
                                                          "T\u00e9l\u00e9charger")))
                )
            )),

        tabPanel(
            title="G\u00e9notyper",
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
                         p("Auteurs: Timoth\u00e9e Flutre, Jacques David"),
                         p("Copyright 2014-2017 INRA, Montpellier SupAgro")
                     ),
                     mainPanel(
                         p("L'acc\u00e8s au g\u00e9nome permet de produire de nombreuses donn\u00e9es de g\u00e9notypage sur un nombre important d'individus qui peuvent \u00eatre ph\u00e9notyp\u00e9s \u00e0 des caract\u00e8res d'int\u00e9r\u00eat agronomique. L'association entre marqueurs et caract\u00e8res fait souvent l'impasse sur les g\u00e8nes ayant de faibles effets, ceux-ci pouvant pourtant d\u00e9terminer une part importante de la variation d'origine g\u00e9n\u00e9tique, et donc du  progr\u00e8s g\u00e9n\u00e9tique envisageable. La pr\u00e9diction g\u00e9nomique des valeurs des plantes, puis l'utilisation de ces pr\u00e9dictions dans un programme de s\u00e9lection, est donc un des enjeux majeurs."),
                         br(),
                         p("Il n'est pas facile d'enseigner ces nouveaux aspects car, d'une part, d'importants pr\u00e9-requis sont n\u00e9cessaires, et d'autre part, il n'existe pas de m\u00e9thodologie d\u00e9finitive valable pour tous mod\u00e8les biologiques, avec des enjeux variables en terme de programmation de progr\u00e8s g\u00e9n\u00e9tique et de sortie vari\u00e9tale. C'est un domaine tr\u00e8s actif de recherche avec de nombreuses innovations en statistique et g\u00e9n\u00e9tique."),
                         p(" Notre ambition est de confronter les \u00e9tudiants \u00e0 certains aspects th\u00e9oriques et pratiques dans une dynamique d'apprentissage actif bas\u00e9e sur la prise en main pratique d'outils et sur l'utilisation d'un \"jeu s\u00e9rieux\"."),
                         br(),
                         p("Les documents p\u00e9dagogiques sont disponibles sur ",
                           a("GitHub.",
                             href="https://github.com/timflutre/atelier-prediction-genomique"))
                     )
                 ))
    )
)
