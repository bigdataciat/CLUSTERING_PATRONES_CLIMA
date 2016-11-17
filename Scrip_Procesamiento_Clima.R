
################################################################################
###-------------------Cargar paquetes, rutas y formatos--------------------#####


libs= c("reshape2","dtwclust","dtw","gtools","agricolae")

lapply(libs, require, character.only=T)

#Especificar ubicación de carpeta de trabajo y leer scripts de funciones

setwd("D:/GIT_HUB_REPOSITORIOS/CLUSTERING_PATRONES_CLIMA")

source('funciones_cluster_temporal.R')

#Lectura de bases de datos con eventos de producción.

base_eventos <- read.csv("DATOS/maiz_cluster_cerete.csv",row.names = 1)

base_eventos$FECHA_SIEMBRA <- as.Date(base_eventos$FECHA_SIEMBRA,"%m/%d/%Y")

base_eventos$FECHA_COSECHA <- as.Date(base_eventos$FECHA_COSECHA,"%m/%d/%Y")
    
#Lectura de datos climaticas. Dos posibilidad de formatos. (1 Separado 2 Unidos)

#Separados

list.files("DATOS/DIVIDIDOS")

nombresClima <- c("ESOL","RAIN","RHUM","TMAX","TMIN")

listDatosClimaticos <-
  lapply(list.files("DATOS/DIVIDIDOS",full.names = T),read.table,header=T)

lapply(listDatosClimaticos,head)

DatosClimaticos <- unifDatos(listDatosClimaticos,nombresClima)

DatosClimaticos$DATE <- as.Date(DatosClimaticos$DATE)

#Unidos

DatosClimaticos <- read.csv("DATOS/13075030_JOINT_CERETE.csv") 

#Acomodar formato de fecha

head(DatosClimaticos)

DatosClimaticos$DATE <- as.Date(DatosClimaticos$DATE,"%m/%d/%Y")

################################################################################
#########------------Division, listas y estandarización ------------############


#REVISAR EL ID QUE SE ESTÁ PERDIENDO


ClimaEventos <-
joinEventsClim(climStat = DatosClimaticos,cropData= base_eventos,datCS = "DATE",
               sowDat = "FECHA_SIEMBRA",harvDat = "FECHA_COSECHA")

head(ClimaEventos)

#Transformación de las variables con eventos climáticos

ClimaEventosTransf <- ClimaEventos

ClimaEventosTransf$LOG10RAIN <- ClimaEventosTransf$RAIN 

ClimaEventosTransf$LOG10RAIN[ClimaEventosTransf$LOG10RAIN==0] <- 0.05

ClimaEventosTransf$LOG10RAIN <- log10(ClimaEventosTransf$LOG10RAIN)

varsI <- c("ESOL","LOG10RAIN","RHUM", "TMAX","TMIN")

#Listas de eventos climaticos normalizados

procesData(climEvent = ClimaEventosTransf,idVar = "EVENT",NormMethod = 2,
           vars = varsI)

load("listClimatEvent.RData")

#Se convierte en series de tiempo multivariadas

tsnleventsN <- lapply(evenN,ts)


################################################################################
#######-----------------------ANALISIS CLUSTER---------------------------#######

distAllMatrix <- distDtwMV(tsnleventsN)

hClustEvents <- hirarCluster(distAllMatrix)

save(distAllMatrix,file = "distMatrixCluster.RData")

IdEvent <- names(evenF)

eventosClasificados <- data.frame(base_eventos,Cluster=hClustEvents)

write.csv(eventosClasificados,"eventsClasificated.csv")


################################################################################
##########------------Generar descriptivos de los cluster----------#############

rsmnClustClima(evenF,hClustEvents,varsI)

m=ggplot(eventosClasificados,aes(FECHA_COSECHA,Rendimiento, 
  colour=as.factor(Cluster)))+geom_point()+theme_bw()+xlab("Fecha Cosecha")+
    ylab("Rendimiendo")+guides(colour=guide_legend(title="Cluster"))

ggsave("RendimientoFechaCluster.png",m,width =6 ,height = 3)



h <- krusk.boxplot(baseConComparacion=eventosClasificados,vary = "Rendimiento",
              varx = "Cluster",
              ylabs=expression(paste("Maize Yield (kg.",ha^-1,")")),maxVar=0)

png("boxplot_numDatosRend.png",height = 350, width = 600)
print(h)
dev.off()

