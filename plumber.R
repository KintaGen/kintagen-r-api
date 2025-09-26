# plumber.R
library(plumber)

#* @apiTitle My First R API
#* @apiDescription An API to demonstrate running R on DigitalOcean.

#* Echo back the message
#* @param msg The message to echo
#* @get /echo
function(msg="") {
  list(msg = paste0("The server received this message: '", msg, "'"))
}

#* Just say hello
#* @get /hello
function() {
  list(message = "Hello from your R API hosted on DigitalOcean!")
}