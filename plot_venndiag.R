library(VennDiagram)


grid.newpage()
draw.pairwise.venn(2067, 1778, 1761, category = c("Two mixure", "TDA"),
                   lty = rep("blank", 2),
                   fill = c("light blue", "pink"),
                   alpha = rep(0.5, 2),
                   #cat.pos = c(0,0),
                   #cat.dist = rep(0.025, 2)
                   )

grid.newpage()
draw.pairwise.venn(575, 487, 484, category = c("Two mixure", "TDA"),
                   lty = rep("blank", 2),
                   fill = c("light blue", "pink"),
                   alpha = rep(0.5, 2),
                   #cat.pos = c(0,0),
                   #cat.dist = rep(0.025, 2)
                   )


