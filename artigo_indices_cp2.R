#### INDICES (sieeenta!)

# exclusão das parcelas que contem somente um ou dois individuo - para calcular os indices é necessario no mínimo 3 indivíduos
parc<-data.frame(table(atrib$quad))
parc2<-parc[which(parc$Freq>2),]
atrib2<-atrib[atrib$quad %in% parc2$Var1,]

# selecionando as subparcelas coletadas na planilha geral do peic
atrib.quad<-sort(unique(atrib2$quad))
censo2<-censo[censo$quad %in% atrib.quad,]

## Calculando a abundancia das especies por subparcela
abund<-data.frame(table(censo2$species, censo2$quad))
parcela<-sort(unique(censo2$quad))
for(i in 1:length(parcela))
{abund$rel[abund$Var2==parcela[i]]<-(abund$Freq[abund$Var2==parcela[i]]/
                                       sum(abund$Freq[abund$Var2==parcela[i]]))*100}

abund$Var2<-as.character(abund$Var2)
abund$Var1<-as.character(abund$Var1)

# Selecionando as subparcelas que possuem 80% de abundância das espécies coletadas
corte<-data.frame(parcela, m.spp=NA)
especies<-as.factor(unique(atrib2$species))

for(i in 1:length(parcela))
{corte$m.spp[corte$parcela==parcela[i]]<-
   sum(abund$rel[abund$Var2==parcela[i] & abund$Var1 %in% especies])}

parcelas.c<-corte$parcela[which(corte$m.spp>80)]

# separando as parcelas selecionadas pelo corte de 80% de abundancia
atrib.corte<-atrib2[atrib2$quad %in% parcelas.c,]

## porcentagem de parcelas em cada tipo de solo
quad1<-data.frame(quad=unique(atrib.corte$quad))
quad1$soil<-atrib.corte$soil[match(quad1$quad,atrib.corte$quad)]
table(quad1$soil) 

# calculando os indices
parcela2<-unique(atrib.corte$quad)
indices<-list()
#require("FD")
for(i in 1:length(parcela2))
{
  atrib.m<-atrib.corte[atrib.corte$quad==parcela2[i],c(29,32,34:36,7)]
  atrib.m2<-aggregate(atrib.m[,1:5],by=list(atrib.m$species),mean)
  atrib.m3<-as.matrix(atrib.m2[,2:6])
  rownames(atrib.m3)<-atrib.m2$Group.1
  
  especies2<-rownames(atrib.m3)
  abund.m<-as.matrix(t(abund$rel[abund$Var2==parcela2[i] & abund$Var1 %in% especies2]))
  colnames(abund.m)<-especies2
  indices[[i]]<-dbFD(atrib.m3, abund.m)
}
head(indices)

###############################################################
##extraindo os valores dos indices do objeto list
names(indices)<-parcela2
indices2<-data.frame( fric=NA, feve=NA, fdiv=NA, cwm.af=NA, cwm.lt=NA,
                      cwm.ldmc=NA, cwm.dm=NA, cwm.sla=NA, quad=parcela2)

indices2$soil<-atrib.corte$soil[match(indices2$quad,atrib.corte$quad)]

for(i in 1:length(parcela2))
{
  indices2$fric[indices2$quad==parcela2[i]]<-indices[[i]][3]
  indices2$feve[indices2$quad==parcela2[i]]<-indices[[i]][5]
  indices2$fdiv[indices2$quad==parcela2[i]]<-indices[[i]][6]
  indices2$cwm.af[indices2$quad==parcela2[i]]<-indices[[i]][[9]][1]
  indices2$cwm.lt[indices2$quad==parcela2[i]]<-indices[[i]][[9]][2]
  indices2$cwm.ldmc[indices2$quad==parcela2[i]]<-indices[[i]][[9]][3]
  indices2$cwm.dm[indices2$quad==parcela2[i]]<-indices[[i]][[9]][4]
  indices2$cwm.sla[indices2$quad==parcela2[i]]<-indices[[i]][[9]][5]
}
head(indices2)

## transformando a classe list em numeric
indices3<- matrix(NA, nrow = dim(indices2)[1], ncol = 8)
for (i in 1:8){ indices3[,i] <- c(as.numeric(indices2[[i]]))}
mode(indices3)
rownames(indices3)<-indices2$soil
colnames(indices3)<-c("fric", "feve", "fdiv", "cwm.af", "cwm.lt",
                      "cwm.ldmc","cwm.dm", "cwm.sla")

###################
###### PCA ########
###################
sum(is.na(indices3))
traits.scaled<-apply(indices3, MARGIN=2, scale)
pc=princomp(traits.scaled, na.rm=T)
biplot(pc)
summary(pc)

# separando os scores para fazer graficos e analises
pc.scores<-pc$scores
pc.scores2<-data.frame(pc.scores)
pc.scores2$soil<-rownames(indices3)
colnames(pc.scores2)<-colnames(indices2)[c(1:8,10)]
head(pc.scores2)

##grafico marcel
pc.loadings<-pc$loadings
pc.var<-matrix((pc$sdev^2/sum(pc$sdev^2))*100,nrow=1)
colnames(pc.var)<-c("cp1","cp2","cp3","cp4","cp5","cp6","cp7", "cp8")

pdf("pca_indices.pdf")
for(i in 1:8)
{for(j in 1:8)
{fbiplot(data=pc.scores2,  x=i , y=j, group= 9,lds=pc.loadings,vs=pc.var,
        arr=3.5,txt=1.2, tam=8)
}}
dev.off()

#####################################
### Selecao de modelos ##############
#####################################
lm01<-lm(Comp.1~soil, data=pc.scores2)
lm02<-update(lm01, .~. -soil)
anova(lm01, lm02) # solo é significativo
AIC(lm01,lm02)
qqnorm(pc.scores2[,1])
qqline(pc.scores2[,1])

lm03<-lm(Comp.2~soil, data=pc.scores2)
lm04<-update(lm03, .~. -soil)
anova(lm03, lm04) # solo não é significativo
qqnorm(pc.scores2[,2])
qqline(pc.scores2[,2])

lm05<-lm(Comp.3~soil, data=pc.scores2)
lm06<-update(lm05, .~. -soil)
anova(lm05, lm06) # solo não é significativo
qqnorm(pc.scores2[,3])
qqline(pc.scores2[,3])

lm07<-lm(Comp.4~soil, data=pc.scores2)
lm08<-update(lm07, .~. -soil)
anova(lm07, lm08) # solo não é significativo
qqnorm(pc.scores2[,4])
qqline(pc.scores2[,4])

lm09<-lm(Comp.5~soil, data=pc.scores2)
lm10<-update(lm09, .~. -soil)
anova(lm09, lm10) # solo não é significativo
qqnorm(pc.scores2[,5])
qqline(pc.scores2[,5])

lm11<-lm(Comp.6~soil, data=pc.scores2)
lm12<-update(lm11, .~. -soil)
anova(lm11, lm12) # solo não é significativo
qqnorm(pc.scores2[,6])
qqline(pc.scores2[,6])

lm13<-lm(Comp.7~soil, data=pc.scores2)
lm14<-update(lm13, .~. -soil)
anova(lm13, lm14) # solo não é significativo
qqnorm(pc.scores2[,7])
qqline(pc.scores2[,7])

lm15<-lm(Comp.8~soil, data=pc.scores2)
lm16<-update(lm15, .~. -soil)
anova(lm15, lm16) # solo não é significativo
qqnorm(pc.scores2[,8])
qqline(pc.scores2[,8])
