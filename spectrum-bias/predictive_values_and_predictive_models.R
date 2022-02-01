n <- 1e6
p1 <- rbeta(n, 1, 4)
y1 <- rbinom(n, 1, p1)
p2 <- rbeta(n, 1/2, 1/2)
y2 <- rbinom(n, 1, p2)

par(mfrow = c(1,2))
hist(p1, 
     xlab = "Prob. of disease", 
     main = "Low signal")
text(
  0.5, 1.5e5,
  paste0(
    "\nThreshold: 69% \u2012 ",
    "PPV: ", round(mean(y1[p1 > .69])*100), "%, ",
    "NPV: ", round(mean(y1[p1 < .69]==0)*100), "%\n"
  )
)
segments(x0=0.69,y0=0,x1=0.69,y1=9e4, col="red", lty = 2, lwd = 2)
hist(p2, 
     xlab = "Prob. of disease", 
     main = "High signal")

text(
  0.5, 110000,
  paste0(
    "Threshold: 35% \u2012 ",
    "PPV: ", round(mean(y2[p2 > .35])*100), "%, ",
    "NPV: ", round(mean(y2[p2 < .35]==0)*100), "%\n"
  )
)

segments(x0=0.35, y0=0, x1=0.35, y1=9e4, col="red", lty=2, lwd = 2)
