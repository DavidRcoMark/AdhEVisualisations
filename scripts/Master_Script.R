#Author: Dave Mark
#Purpose: Plot TRIFM Image Data

library(tidyverse)
library(gganimate)
library(ggpubr)
library(matchmaker)
#library(extrafont)
#font_import() #Plots use arial font - if you don't import them it breaks the plots
library(here)

#Read in Everything----

Csv <- list.files(path = here(), recursive = T, pattern = "csv") #Lists csv files
PhotonList <- lapply(Csv, read.csv) #Reads the files in as a list object

Names <- data.frame("ID"=Csv, "Num"=1:15) #Creates a "Dictionary" so we can name things down the line

DF <- bind_rows(PhotonList, .id="FileName") #Turns list into one table
DF$Name <- match_vec(x = DF$FileName, dictionary = Names, from = "Num", to = "ID") #Links ID to each row in the table


#Add fluorophore and drug conc ----
DF <- DF |> mutate("Fluorophore"=case_when(grepl("Green", DF$Name)~"J549",
                                     grepl("Red", DF$Name)~"J646",
                                     grepl("Yellow", DF$Name)~"Merge")) |>
  mutate("Compound_Ratio"=case_when(grepl("1.1 compound", DF$Name)~"1:1",
                                    grepl("1.10 compound", DF$Name)~"1:10", .default = "1:0"))

DF <- DF |> pivot_longer(cols = 4:16)

DF$name <- gsub("T", "", DF$name)
DF$name <- as.numeric(DF$name)
colnames(DF)[7]="Time(min)"


#Make the videos ----
##S1----
MovieS1 <- DF |> filter(Compound_Ratio=="1:0") |>
  ggplot(aes(x=Bins, y=value, col=`Time(min)`, fill=`Time(min)`))+
  stat_summary(geom="col", position = "identity", alpha=0.7)+
  stat_summary(geom="errorbar")+
  transition_time(`Time(min)`)+
  labs(title = "{frame_time} Minutes", subtitle = "Fluorophore")+
  facet_wrap(~factor(Fluorophore, levels = c("J646", "J549", "Merge")))+
  labs(x = "Maximum intensity (photons)", y = "Number of spirosomes")+
  scale_color_continuous(type = "viridis")+
  scale_fill_continuous(type = "viridis")+
  scale_y_continuous(trans="pseudo_log", breaks = c(0,1e1,1e2,1e3,1e4))+
  theme_bw()+
  theme(legend.position = "none",
        text = element_text(family = "Arial"),
        axis.text = element_text(size = 8),
        strip.background = element_rect(color="transparent", fill="white"),
        strip.text = element_text(size=10),
        plot.title = element_text(size=14),
        plot.subtitle = element_text(size=12, hjust = 0.5))+
  coord_cartesian(y=c(1, 1.1e4), x=c(0, 2.67e3))

histS1 <- animate(MovieS1, renderer = gifski_renderer(), fps = 2, nframes = 13,
                  width=20,height=8, units="cm", res=150)
histS1
##S2----
MovieS2 <- DF |> filter(Fluorophore=="J549") |>
  ggplot(aes(x=Bins, y=value, col=`Time(min)`, fill=`Time(min)`))+
  stat_summary(geom="col", position = "identity", alpha=0.7)+
  stat_summary(geom="errorbar")+
  transition_time(`Time(min)`)+
  labs(title = "{frame_time} Minutes", subtitle = "AdhE:ME0054 (J549)")+
  facet_wrap(~factor(Compound_Ratio, levels = c("1:0", "1:1", "1:10")))+
  labs(x = "Maximum intensity (photons)", y = "Number of spirosomes")+
  scale_color_continuous(type = "viridis")+
  scale_fill_continuous(type = "viridis")+
  scale_y_continuous(trans="pseudo_log", breaks = c(0,1e1,1e2,1e3,1e4))+
  theme_bw()+
  theme(legend.position = "none",
        text = element_text(family = "Arial"),
        axis.text = element_text(size = 8),
        strip.background = element_rect(color="transparent", fill="white"),
        strip.text = element_text(size=10),
        plot.title = element_text(size=14),
        plot.subtitle = element_text(size=12, hjust = 0.5))+
  coord_cartesian(y=c(1, 1.1e4), x=c(0, 2.67e3))

histS2 <- animate(MovieS2, renderer = gifski_renderer(), fps = 2, nframes = 13,
                  width=20,height=8, units="cm", res=150)
histS2

##S3----
MovieS3 <- DF |> filter(Fluorophore=="J646") |>
  ggplot(aes(x=Bins, y=value, col=`Time(min)`, fill=`Time(min)`))+
  stat_summary(geom="col", position = "identity", alpha=0.7)+
  stat_summary(geom="errorbar")+
  transition_time(`Time(min)`)+
  labs(title = "{frame_time} Minutes", subtitle = "AdhE:ME0054 (J646)")+
  facet_wrap(~factor(Compound_Ratio, levels = c("1:0", "1:1", "1:10")))+
  labs(x = "Maximum intensity (photons)", y = "Number of spirosomes")+
  scale_color_continuous(type = "viridis")+
  scale_fill_continuous(type = "viridis")+
  scale_y_continuous(trans="pseudo_log", breaks = c(0,1e1,1e2,1e3,1e4))+
  theme_bw()+
  theme(legend.position = "none",
        text = element_text(family = "Arial"),
        axis.text = element_text(size = 8),
        strip.background = element_rect(color="transparent", fill="white"),
        strip.text = element_text(size=10),
        plot.title = element_text(size=14),
        plot.subtitle = element_text(size=12, hjust = 0.5))+
  coord_cartesian(y=c(1, 1.1e4), x=c(0, 2.67e3))

histS3 <- animate(MovieS3, renderer = gifski_renderer(), fps = 2, nframes = 13,
                  width=20,height=8, units="cm", res=150)
histS3
##S4----
MovieS4 <- DF |> filter(Fluorophore=="Merge") |>
  ggplot(aes(x=Bins, y=value, col=`Time(min)`, fill=`Time(min)`))+
  stat_summary(geom="col", position = "identity", alpha=0.7)+
  stat_summary(geom="errorbar")+
  transition_time(`Time(min)`)+
  labs(title = "{frame_time} Minutes", subtitle = "AdhE:ME0054 (Merge)")+
  facet_wrap(~factor(Compound_Ratio, levels = c("1:0", "1:1", "1:10")))+
  labs(x = "Maximum intensity (photons)", y = "Number of spirosomes")+
  scale_color_continuous(type = "viridis")+
  scale_fill_continuous(type = "viridis")+
  scale_y_continuous(trans="pseudo_log", breaks = c(0,1e1,1e2,1e3,1e4))+
  theme_bw()+
  theme(legend.position = "none",
        text = element_text(family = "Arial"),
        axis.text = element_text(size = 8),
        strip.background = element_rect(color="transparent", fill="white"),
        strip.text = element_text(size=10),
        plot.title = element_text(size=14),
        plot.subtitle = element_text(size=12, hjust = 0.5))+
  coord_cartesian(y=c(1, 1.1e4), x=c(0, 2.67e3))

histS4 <- animate(MovieS4, renderer = gifski_renderer(), fps = 2, nframes = 13,
                  width=20,height=8, units="cm", res=150)
histS4

##Save the animations ----

anim_save("./figures/SupportingMovies/MovieS1.gif", histS1)
anim_save("./figures/SupportingMovies/MovieS2.gif", histS2)
anim_save("./figures/SupportingMovies/MovieS3.gif", histS3)
anim_save("./figures/SupportingMovies/MovieS4.gif", histS4)

#Draw Fig 2D ----

Fig2D <- DF |> filter(Compound_Ratio=="1:0") |> filter(Fluorophore=="Merge") |>
  ggplot(aes(x=Bins, y=value, col=`Time(min)`, fill=`Time(min)`))+
  stat_summary(geom="col", position = "identity", alpha=0.7)+ #Draws columns
  stat_summary(geom="errorbar")+ #Draws error bars
  facet_grid(rows=vars(`Time(min)`),)+ #Draws one plot per time limit
  labs(x = "Maximum intensity (photons)", y = "Number of spirosomes")+ #Axis labels
  scale_color_continuous(type = "viridis")+ #Color
  scale_fill_continuous(type = "viridis")+ #Fill
  scale_y_continuous(trans="pseudo_log", breaks = c(0,1e1,1e2,1e3))+ #Manual tweaks to the axis
  theme_bw()+ #Background theme
  coord_cartesian(y=c(1, 1.1e3), x=c(0, 2.67e3))+ #Manually sets the space that can be drawn in
  theme(text=element_text(family="arial"), axis.title = element_text(size=11.5, color="black"), axis.text = element_text(size=9.5, colour = "black"),
        strip.background = element_rect(fill="white"), legend.position = "none")

ggsave(Fig2D, dpi=300, width=5, height = 10, filename="./figures/plots/Fig2D.png")

#Draw Fig S4A----

S4A <- DF |> filter(Fluorophore=="J646") |>
  ggplot(aes(x=Bins, y=value, col=`Time(min)`, fill=`Time(min)`))+
  stat_summary(geom="col", position = "identity", alpha=0.7)+
  stat_summary(geom="errorbar")+
  #transition_time(`Time(min)`)+
  labs(subtitle = "AdhE:ME0054 (J646)")+
  facet_grid(cols = vars(Compound_Ratio), rows = vars(`Time(min)`))+
  labs(x = "Maximum intensity (photons)", y = "Number of spirosomes")+
  scale_color_continuous(type = "viridis")+
  scale_fill_continuous(type = "viridis")+
  scale_y_continuous(trans="pseudo_log", breaks = c(0,1e1,1e2,1e3,1e4))+
  theme_bw()+
  theme(legend.position = "none",
        text = element_text(family = "Arial"),
        axis.text = element_text(size = 8),
        strip.background.x = element_rect(color="transparent", fill="white"),
        strip.background.y = element_rect(fill="white"),
        strip.text = element_text(size=10),
        plot.title = element_text(size=14),
        plot.subtitle = element_text(size=12, hjust = 0.5))+
  coord_cartesian(y=c(1, 1.1e4), x=c(0, 2.67e3))

ggsave(S4A, filename="./figures/plots/S4A.png", dpi=300,
       width = 17.4, height =21, units = "cm")

#Draw Fig S4B----
S4C <- DF |> filter(Fluorophore=="J549") |>
  ggplot(aes(x=Bins, y=value, col=`Time(min)`, fill=`Time(min)`))+
  stat_summary(geom="col", position = "identity", alpha=0.7)+
  stat_summary(geom="errorbar")+
  #transition_time(`Time(min)`)+
  labs(subtitle = "AdhE:ME0054 (J549)")+
  facet_grid(cols = vars(Compound_Ratio), rows = vars(`Time(min)`))+
  labs(x = "Maximum intensity (photons)", y = "Number of spirosomes")+
  scale_color_continuous(type = "viridis")+
  scale_fill_continuous(type = "viridis")+
  scale_y_continuous(trans="pseudo_log", breaks = c(0,1e1,1e2,1e3,1e4))+
  theme_bw()+
  theme(legend.position = "none",
        text = element_text(family = "Arial"),
        axis.text = element_text(size = 8),
        strip.background.x = element_rect(color="transparent", fill="white"),
        strip.background.y = element_rect(fill="white"),
        strip.text = element_text(size=10),
        plot.title = element_text(size=14),
        plot.subtitle = element_text(size=12, hjust = 0.5))+
  coord_cartesian(y=c(1, 1.1e4), x=c(0, 2.67e3))

ggsave(S4C, filename="./figures/plots/S4B.png", dpi=300,
       width = 17.4, height =21, units = "cm")

#Draw Fig S4C----
S4C <- DF |> filter(Fluorophore=="Merge") |>
  ggplot(aes(x=Bins, y=value, col=`Time(min)`, fill=`Time(min)`))+
  stat_summary(geom="col", position = "identity", alpha=0.7)+
  stat_summary(geom="errorbar")+
  #transition_time(`Time(min)`)+
  labs(subtitle = "AdhE:ME0054 (J549)")+
  facet_grid(cols = vars(Compound_Ratio), rows = vars(`Time(min)`))+
  labs(x = "Maximum intensity (photons)", y = "Number of spirosomes")+
  scale_color_continuous(type = "viridis")+
  scale_fill_continuous(type = "viridis")+
  scale_y_continuous(trans="pseudo_log", breaks = c(0,1e1,1e2,1e3,1e4))+
  theme_bw()+
  theme(legend.position = "none",
        text = element_text(family = "Arial"),
        axis.text = element_text(size = 8),
        strip.background.x = element_rect(color="transparent", fill="white"),
        strip.background.y = element_rect(fill="white"),
        strip.text = element_text(size=10),
        plot.title = element_text(size=14),
        plot.subtitle = element_text(size=12, hjust = 0.5))+
  coord_cartesian(y=c(1, 1.1e4), x=c(0, 2.67e3))

ggsave(S4C, filename="./figures/plots/S4C.png", dpi=300,
       width = 17.4, height =21, units = "cm")
