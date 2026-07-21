# Figure 1: One Health transmission dynamics of RVFV in Africa

getwd()
setwd("~/Desktop/Manuscripts/Review_Article/")

if (!require(DiagrammeR)) install.packages("DiagrammeR")
if (!require(DiagrammeRsvg)) install.packages("DiagrammeRsvg")
if (!require(rsvg)) install.packages("rsvg")

library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

g <- grViz("
digraph rvfv_one_health_clean {

  graph [
    layout = dot,
    rankdir = TB,
    bgcolor = white,
    pad = 0.3,
    nodesep = 0.6,
    ranksep = 0.8,
    splines = curved,
    overlap = false
  ]

  node [
    shape = box,
    style = 'rounded,filled',
    color = black,
    penwidth = 1.2,
    fontname = Helvetica,
    fontsize = 16,
    margin = 0.18
  ]

  edge [
    color = black,
    penwidth = 1.4,
    fontname = Helvetica,
    fontsize = 12,
    arrowsize = 0.8
  ]

  # Top layer
  env [
    label = 'Environmental drivers\\nHeavy rainfall\\nFlooding\\nIrrigation and vegetation\\nClimate variability',
    fillcolor = '#DCEEF8'
  ]

  # Middle layer
  mosquito [
    label = 'Mosquito vectors\\nAedes spp.\\nCulex spp.',
    fillcolor = '#F7E39C'
  ]

  virus [
    label = 'RVFV',
    shape = ellipse,
    style = filled,
    fillcolor = '#D9DDE3',
    fontsize = 20,
    width = 1.5,
    height = 0.7
  ]

  livestock [
    label = 'Livestock hosts\\nSheep\\nGoats\\nCattle',
    fillcolor = '#CFEFD6'
  ]

  humans [
    label = 'Humans\\nMosquito bites\\nContact with blood, tissues,\\nand aborted materials',
    fillcolor = '#F6D1D1'
  ]

  # Bottom layer
  onehealth [
    label = 'One Health response\\nDiagnostics   |   Genomic surveillance   |   Vaccination   |   Data sharing',
    fillcolor = '#E6DDF8',
    width = 7.5
  ]

  # Rank control
  { rank = same; mosquito; virus; livestock }
  { rank = same; humans }
  { rank = same; onehealth }

  # Main biological flow
  env -> mosquito [
    label = 'increases vector abundance',
    color = '#2563EB',
    fontcolor = '#2563EB'
  ]

  mosquito -> livestock [
    label = 'vector-borne transmission'
  ]

  mosquito -> humans [
    label = 'mosquito bite'
  ]

  livestock -> humans [
    label = 'occupational exposure'
  ]

  # Virus anchor
  mosquito -> virus [
    style = dashed,
    arrowhead = none,
    color = gray60
  ]

  virus -> livestock [
    style = dashed,
    arrowhead = none,
    color = gray60
  ]

  # One Health support links
  onehealth -> mosquito [
    style = dotted,
    color = '#7C3AED',
    arrowhead = none
  ]

  onehealth -> livestock [
    style = dotted,
    color = '#7C3AED',
    arrowhead = none
  ]

  onehealth -> humans [
    style = dotted,
    color = '#7C3AED',
    arrowhead = none
  ]
}
")

# View figure in RStudio
g

# Export as high-resolution PNG to Desktop
svg_txt <- export_svg(g)

rsvg_png(
  charToRaw(svg_txt),
  file = '~/Desktop/Manuscripts/Review_Article/Figure1_RVFV_OneHealth_clean.png',
  width = 2200,
  height = 1600
)
