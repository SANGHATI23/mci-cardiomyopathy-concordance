# Shiny Deployment Notes

Recommended deployment target:

- shinyapps.io for manuscript review
- persistent URL before manuscript submission

Before deployment:

1. Confirm results/resource_tables/MCI_BROWSER_READY_DISPLAY_TABLE_v0_2.csv exists.
2. Confirm shiny_app/app.R loads the table correctly.
3. Run locally with shiny::runApp('shiny_app').
4. Deploy using rsconnect.

Example R commands:

install.packages('rsconnect')
library(rsconnect)

rsconnect::setAccountInfo(
  name = 'YOUR_SHINYAPPS_NAME',
  token = 'YOUR_TOKEN',
  secret = 'YOUR_SECRET'
)

rsconnect::deployApp('shiny_app')

Do not commit tokens or secrets to GitHub.
