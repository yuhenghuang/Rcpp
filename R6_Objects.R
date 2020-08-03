library(Rcpp)
library(R6)

# optional
# redirect cache directory of rcpp
options(rcpp.cache.dir = "D:\\Programming\\Advanced_R\\Rcpp_temp")

sourceCpp("R6_Objects.cpp")


#' R6 class that defines a person
Person <- R6::R6Class("Person",
                      public = list(
                        #' @description 
                        #' Constructor of the 'Person' class
                        #' @param name a string that defines this person's name
                        #' @param id an integer that defines this person's id
                        #' @return A new 'Person' object
                        initialize = function(name, id){
                          private$name <- name
                          private$id <- id
                        },
                        
                        get_name = function(){return(private$name)},
                        
                        get_id = function(){return(private$id)},
                        
                        #' @description 
                        #' Gives an item to the person
                        #' @param item a string that defines the item
                        give_item = function(item){private$item <- item},
                        
                        #' @description 
                        #' A public function that calls a private one
                        good_foo = function(){
                          return(paste0("Wrapped inside: {", private$bad_foo(), "}"))
                        }
                      ),
                      private = list(
                        #' @field name the name of the person
                        name = NULL,
                        #' @field id the id of the person
                        id = NULL,
                        #' @field item some item that the person has
                        item = NULL,
                        
                        #' @description 
                        #' A private function that should not be called from outside the object
                        bad_foo = function(){return("This is a private function")}
                      )
)


names <- c("Jake", "Anne", "Alex")
ids <- 1:3

res <- initialize_list(names, ids)

print(res)


tryCatch({res <- access_method("bad_foo")},
         error=function(cond){print(paste0("Exception raised: ", cond))})

res <- access_method("good_foo")

print(res)



p <- Person$new("Jake", 1)

give_item(p, "shovel")

print(p)
