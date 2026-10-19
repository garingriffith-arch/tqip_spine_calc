# SPINE-TRACT Shiny entry point
# TBI-TRACT-style interface, SPINE-TRACT locked 2020–2024 deployment artifacts.
# The application implementation and inference logic are in app_core.R.
suppressPackageStartupMessages(library(shiny))
app <- source("app_core.R", local = TRUE)$value
app
