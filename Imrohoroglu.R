library(Rcpp)
library(ggplot2)

sourceCpp("Imrohoroglu.cpp")

st <- Sys.time()

res <- imrohoroglu(0.995, 1.5, 0.25)

total <- Sys.time() - st
print(total)

ggplot(mapping = aes(x = res$a_grid)) +
  geom_line(mapping = aes(y=res$demo[,1]), color = "blue") +
  geom_line(mapping = aes(y=res$demo[,2]), color = "red")


ggplot(mapping = aes(x = res$a_grid)) +
  geom_line(mapping = aes(y=res$demo[,3]), color = "blue") +
  geom_line(mapping = aes(y=res$demo[,4]), color = "red")


ggplot(mapping = aes(x = res$a_grid)) +
  geom_line(mapping = aes(y=res$policy[,1]/37.5), color = "blue") +
  geom_line(mapping = aes(y=res$policy[,2]/37.5), color = "red") +
  geom_abline(slope = 1, color="black")


ggplot(mapping = aes(x = res$a_grid)) +
  geom_line(mapping = aes(y=res$policy[,3]/37.5), color = "blue") +
  geom_line(mapping = aes(y=res$policy[,4]/37.5), color = "red") +
  geom_abline(slope = 1, color="black")
