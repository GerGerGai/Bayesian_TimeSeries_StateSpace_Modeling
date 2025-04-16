getwd()
setwd("RScripts")

data <- read.table("../earthquake_california/earthquake.txt",
                   header = FALSE,
                   skip = 2,
                   sep = "",
                   stringsAsFactors = FALSE,
                   fill = TRUE)


colnames(data) <- c("Date", "Time", "Lat", "Lon", "Depth", "Mag", "Magt",
                    "Nst", "Gap", "Clo", "RMS", "SRC", "Event_ID")


data <- na.omit(data)

write.csv(data, file = "../cleaned_earthquake_data.csv", row.names = FALSE)
