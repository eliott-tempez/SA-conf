#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# To modify
#-------------------------------------------------------------------------------

folder_out = "saconf_out"
seuilVar = 1.5


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Variables
#-------------------------------------------------------------------------------

source(paste(folder_out, "params.R", sep="/"))

#-----------------------------------------
# File input
#-----------------------------------------
Listpdb = as.character(read.table(file1)[,1])

fileLS = paste(pathsaconf,"SL_alignment.fasta2_wo_header.aln",sep="/")
LS = read.table(fileLS)

fileAA = paste(pathsaconf,"AA_alignment.fasta2_wo_header.aln",sep="/")
AA = read.table(fileAA)

file.correp = paste(pathsaconf, "corresp_Swss_align.csv",sep="/")
fileC = read.table(file.correp,sep=",",header=T)

#-----------------------------------------
# File output
#-----------------------------------------
file.AAalign = paste(pathFig, "AA-alignment.pdf", sep="/")
file.SLalign = paste(pathFig, "SL-alignment.pdf", sep="/")
file.neq = paste(pathFig, "Neq_graph.pdf", sep="/")
file.out = paste(pathsaconf,"Position_description.csv",sep="/")
file_out1 =paste(pathsaconf, "Count_position_type.txt",sep="/")
filePos_out = paste(pathsaconf, "Correspondence_positions.csv",sep="/")
file.mut = paste(pathsaconf, "Mutation_res.txt",sep="/")
file.SLVar = paste(pathsaconf, "Structural_Variable_position_res.txt",sep="/")



#-----------------------------------------
# variable definition
#-----------------------------------------
#AA
AlphabetAA = c("-","A","I","L","M","F","W","V","R","K","N","Q","S","T","C","E","D","G","H","Y","P")
VectcolAA = c("purple4","blue4","skyblue4","blue1","skyblue","paleturquoise2","lightblue1","red1","red3",
              "olivedrab1","green","green3","olivedrab4","pink","magenta1","magenta3","orange","cyan1",
              "cyan4","yellow2")
names(VectcolAA) =  AlphabetAA[2:length(AlphabetAA)]      


#SL
Alphabet=c("-","a","A","V","W","Z","B","C","D","E","O","S","R","Q","I","F","U","P","H","G","Y","J","K",
           "L","N","M","T","X")
VectcolSL=c("#FF0000FF", "#FF3900FF", "#FF7100FF" ,"#FFAA00FF", "#FFE300FF" , "#FF0071FF","#FF0039FF" ,"#FF00AAFF" ,   "#FF00E3FF" , "#E300FFFF", "#AA00FFFF", "#7100FFFF", "#3900FFFF", "#0000FFFF", "#0039FFFF" , "#0071FFFF","#00AAFFFF","#00E3FFFF" ,"#00FFE3FF", "#00FFAAFF","#00FF71FF", "#00FF39FF","#39FF00FF","#00FF00FF" ,"#71FF00FF", "#AAFF00FF","#E3FF00FF" )
VectcolSL=c("red4","red","indianred3", "tomato2","orange","yellow","yellow3","seashell2","mistyrose1","pink","violet","pink3","purple","purple4","lavenderblush4","blue","steelblue","cyan","paleturquoise","powderblue","palegoldenrod",'bisque4',"darkcyan","green","yellowgreen","olivedrab4","khaki4")
names(VectcolSL) =  Alphabet[2:length(Alphabet)]  


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------

f1 = function(vect){
   v1 = table(vect)
   v2 = length(which(names(v1)!="NA"))
   return(v2)
}

f2 = function(vect){
   vect[which(vect=="NA")]="-"
   v1 = paste(names(table(vect)), collapse="")
   v2 = gsub('[AaVW]',"h",v1 )
   v3 = gsub('[LMNTX]',"s",v2 )
   v4 = gsub('[BCDEFGHIJKOPQRSUYZ]',"l",v3)
   v5 = unlist(strsplit(v4,""))
   v5[which(v5=="-")]="NA"
   nbr = f1(v5)
   return(nbr)
}



f21 = function(vect){
   vect[which(vect=="NA")]="-"
   v1 = paste(names(table(vect)), collapse="")
   v2 = gsub('[AaVW]',"h",v1 )
   v3 = gsub('[LMNTX]',"s",v2 )
   v4 = gsub('[BCDEFGHIJKOPQRSUYZ]',"l",v3)
   v5 = unlist(strsplit(v4,""))
   v5[which(v5=="-")]="NA"
   return(v5)
}


f.defSS = function(vect){
   pourc = length(which(vect!="NA")) / length(vect)
   vect = vect[which(vect!="NA")]
   if (pourc >= 0.5){
       nbr = f2(vect)
       listSS = f21(vect)
       if (nbr == 1){
          ss = unique(listSS)
       }else if (nbr == 0){
          ss = "NA"
       }else {
          ss = "v"
       }
    }else{
        ss="NA"
    }
   return(ss)
}



sumFunc = function(vect){
   return(sum(as.numeric(vect), na.rm=T))
}



f.defColSS = function(mm){
    ssCol2 = rep("black", length = dim(mm)[2])
    ssVect = apply(matLS, 2, f.defSS)
    ssCol2[which(ssVect=="v")]="purple"
    ssCol2[which(ssVect=="l")]="gray"
    ssCol2[which(ssVect=="s")]="magenta"
    ssCol2[which(ssVect=="h")]="red"
    ssCol2[which(ssVect=="-")]="white"
    return(ssCol2)
}



computNeq = function(mat){
    neqVect = NULL
    for ( i in 1:dim(mat)[2]){
        vect = mat[which(mat[,i]!="NA"),i]
        pourc = length(vect)/dim(mat)[1]
        if (pourc < 0.50){
              neq = 0
        }else{
             freq = table(vect) / sum(table(vect))
             H = - sum(freq * log(freq))
             neq = exp(H)
        }
        neqVect = c(neqVect, neq)
    }
    return(neqVect)
}




f.nbrGap = function(mat){
    gapVect = NULL
    for ( i in 1:dim(mat)[2]){
        vect = mat[which(mat[,i]=="NA"),i]
        nbr.gap = round(length(vect)/dim(mat)[1],2)
        gapVect = c(gapVect, nbr.gap)
    }
    return(gapVect)
}


f.varSS = function(vect){
   pourc = length(which(vect!="NA")) / length(vect)
   vect = vect[which(vect!="NA")]
   if (pourc >= 0.5){
       nbr = f2(vect)
       listSS = f21(vect)
       if (nbr == 1){
          ltw.v = names(table(listSS))
       }else if (nbr == 0){
          ltw.v = 0     
       }else if (nbr == 3){
          ssU = table(listSS)
          ltw.v = "l-h-s"
       }else if (nbr == 2){
       	  ssU = table(listSS)
          if (( names(ssU)[1] == "h") & (names(ssU)[2] == "l")){
          	ltw.v = "h-l"
          }          
          if (( names(ssU)[1] == "l") & (names(ssU)[2] == "s")){
          	ltw.v = "s-l"
          }          
          if  (( names(ssU)[1] == "h") & (names(ssU)[2] == "s")){
          	ltw.v = "h-s"
          }                                    
       }
    }else{
        ltw.v="-"
    }
   return(ltw.v)
}

		



f.transfoNeqSL = function(neqLSVect, type){
    #type = c("Conser", "WeakVar", "StrongVar")
    if (is.element(type, c("Conser", "WeakVar", "StrongVar")) == FALSE){
       print ("type is not correct")
    }
    if (type == "Conser"){
        seqCompte = rep("0",length=length(neqLSVect))
        seqCompte [which(neqLSVect==1)] = 1
    }

    if (type == "WeakVar"){
        seqCompte = rep("0",length=length(neqLSVect))
        seqCompte [which(neqLSVect < seuilVar)] = 1
        seqCompte [which(neqLSVect <= 1)] = 0
    }

    if (type == "StrongVar"){
        seqCompte = rep("0",length=length(neqLSVect))
        seqCompte [which(neqLSVect>=seuilVar)] = 1
    }

    return(seqCompte)
}




f.ConservAlongSeq = function(seqCompte){
    i = 1
    vtaille = NULL
    vstart = NULL
    while(i <= length(seqCompte)){
        val = seqCompte[i]
        if (val == "1"){
            taille  = 1
            vstart = append(vstart,i)
            while(seqCompte[i+1] == "1"){
                if (i != length(seqCompte)){
                    taille = taille +1
                    i = i+1
                }
                if (i == length(seqCompte)){ 
                   break
                }               
            }
            vtaille = append(vtaille, taille)
        }
        i = i+1
    }
    if (i == length(seqCompte)){
        val = seqCompte[i]
        if (val == "1"){
            taille  = 1
            vtaille = append(vtaille, taille)
        }
    }    
    return(list(vtaille, vstart))
}




f.WriteConservAlongSeq = function(neqLSVect, type="WeakVar", file = file_out1){

    if (is.element(type, c("Conser", "WeakVar", "StrongVar")) == TRUE){
       phrase = c("SL structural conserved",
                  "SL weakly structural variable",
                  "SL strongly structural variable")
       names(phrase) = c("Conser", "WeakVar", "StrongVar")
    
   
        seqCompte = f.transfoNeqSL(neqLSVect, type=type)
        nbrpos.tmp = length(which(seqCompte=="1"))
        if (nbrpos.tmp==0){
            ph0 = paste("\n* No", phrase[type], "in the MSLA", sep=" ")
            write(ph0, file = file_out1, append=T)
        }else{
                ph0 = paste("\n*", nbrpos.tmp,phrase[type], "in the MSLA", sep=" ")
                write(ph0, file = file_out1, append=T)
                l.AlongSeq = f.ConservAlongSeq(seqCompte)
                resultAlongSeq = table(l.AlongSeq[[1]])


                if (is.element("1", names(resultAlongSeq))){
                    phi = paste("    -", resultAlongSeq["1"], phrase[type],"positions are isolated",sep=" ")
                    write(phi, file = file_out1, append=T)    
                    #print(phi)
                    index.startPos = which(l.AlongSeq[[1]]==1)
                    phj = paste("        start positions:", paste(l.AlongSeq[[2]][index.startPos], collapse = " ; "),sep=" ")
                    write(phj, file = file_out1, append=T)    

                }else{
                    phi = paste("    -", 0 , phrase[type],"positions are isolated",sep=" ")
                    write(phi, file = file_out1, append=T)
                    #print(phi)
                }


                for (i in 2:length(resultAlongSeq)){
                    val = resultAlongSeq[i]
                    phi = paste("    -", val, phrase[type], "positions are grouped in regions of ", names(resultAlongSeq)[i], "positions",sep=" ")
                    write(phi, file = file_out1, append=T)
                    #print (phi)

                    index.startPos = which(l.AlongSeq[[1]]==as.numeric(names(resultAlongSeq)[i]))
                    phj = paste("        start positions:", paste(l.AlongSeq[[2]][index.startPos], collapse = " ; "),sep=" ")
                    write(phj, file = file_out1, append=T) 
                }
                
       }
    }else{
       print("type is not correct")
   }
}



f.ExtractPDBRef.Pos = function(ali.pos, folder_out){
   ali =  readLines(paste(folder_out,'position_alignment.fasta2',sep="/"))
   listePos = unlist(strsplit(ali[2],'\\.'))
   return(listePos[ali.pos])
}


f.ExtractPDBRef = function(folder_out){
   ali =  readLines(paste(folder_out,'position_alignment.fasta2',sep="/"))
   protRef = unlist(strsplit(ali[1],'>'))[2]
   return(protRef)
}

f.ExtractPDB.Pos = function(folder_out){
   matPos = NULL
   listProt = NULL
   ali =  readLines(paste(folder_out,'position_alignment.fasta2',sep="/"))
   for (i in 1:length(ali)){
       listetmp = unlist(strsplit(ali[i],''))
       if (is.element(">", listetmp)){
          prot = unlist(strsplit(ali[i],">"))[2]
          listProt = append(listProt, prot)
       }else{
           listePos = unlist(strsplit(ali[i],'\\.'))
           matPos = cbind(matPos, listePos)
       } 

       
   }
   return(list(matPos = matPos, listProt = listProt))
}












#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------


#-----------------------------------------
# Step 1 : compute the number of AA 
# by position
#-----------------------------------------
tailleAA = length(unlist(strsplit(as.character(AA[1,]), "")))
matAA = matrix(NA,ncol = tailleAA, nrow=dim(AA)[1])

for (i in 1:dim(AA)[1]){
   ligne = as.character(AA[i,])
   vect = unlist(strsplit(ligne,""))
   matAA[i,] = vect
}

for (i in 1:dim(matAA)[1]){
   indSel = which(matAA[i,]=="-")
   matAA[i, indSel] = "NA"
}

#recherche le nombre de LS par position
nbrAAbyPos = apply(matAA, 2, f1)

vectAANum = rep(100, length=length(nbrAAbyPos))
vectAANum[which(nbrAAbyPos>1)] = 400

vectAALab = rep(NA, length=length(nbrAAbyPos))
vectAALab[which(nbrAAbyPos>1)] = nbrAAbyPos[which(nbrAAbyPos>1)]


#-----------------------------------------
# Step 2: AA-maxtrix
# visualization of the AA-alignment
#-----------------------------------------

Mat4 = matrix(NA, ncol = ncol(matAA), nrow= nrow(matAA))

for (Ibs in 2:length(AlphabetAA)){
    aa=AlphabetAA[Ibs]
    Mat4[which(matAA==aa)]=as.numeric(Ibs)
}
rownames(Mat4) = Listpdb


pdf(file.AAalign,width=12,paper="USr")
par(mar = c(3, 3.5,1 ,0))   ###c(bottom, left, top, right). Défaut : 5.1 4.1 4.1 2.1
nf = layout(matrix(c(1,2),1,2,byrow=TRUE), c(22.7,4), c(20,20), TRUE)
image(t(Mat4)[,nrow(Mat4):1],ylab="",cex.lab=1.5,cex.main=1.5,col=VectcolAA,axes=F,main="Amino acid alignment")
      #breaks=seq(min(Mat4,na.rm=T),(max(Mat4,na.rm=T)+1))
axis(2,seq(1,0,length=dim(Mat4)[1]),rownames(Mat4),cex.axis=0.6,las=1)
axis(1,seq(0,(dim(Mat4)[2]-1), by=10)/dim(Mat4)[2],seq(1,dim(Mat4)[2], by=10),cex.axis=0.6,las=1)
#axis(4,seq(0,(dim(Mat4)[1]-1), by=5)/dim(Mat4)[1],sort(seq(1,dim(Mat4)[1], by=5),decreasing=T),cex.axis=0.6,las=1)

box()
numero=as.numeric(names(table(Mat4)))
plot(rep(0,length=length(numero)),numero,type="n",axes=F,ylab="",xlab="")
text(rep(0,length=length(numero)),numero,label=AlphabetAA[numero],col=VectcolAA)	
dev.off()



nbgraph = round(dim(Mat4)[1]/50,0)

if (nbgraph>1){
    tailleLim = round(seq(0,dim(Mat4)[1], length=(nbgraph+1)),0)
    for(i in 1:(length(tailleLim)-1)){
        file.AAaligni = paste(pathFig, paste("AA-alignment_", i, ".pdf", sep=""), sep="/")
        pdf(file.AAaligni,width=12,paper="USr")
        par(mar = c(3, 3.5,3 ,0))  
        ssMat = Mat4[(tailleLim[i]+1):tailleLim[(i+1)], ]
        nf = layout(matrix(c(1,2),1,2,byrow=TRUE), c(22.7,4), c(20,20), TRUE)
        image(t(ssMat)[,nrow(ssMat):1],ylab="",cex.lab=1.5,cex.main=1.5,col=VectcolAA,axes=F,main="Amino acid alignment")
        axis(2,seq(1,0,length=dim(ssMat)[1]),rownames(ssMat),cex.axis=0.6,las=1)
        axis(1,seq(0,(dim(ssMat)[2]-1), by=10)/dim(ssMat)[2],seq(1,dim(ssMat)[2], by=10),cex.axis=0.6,las=1)

        box()
        numero=as.numeric(names(table(Mat4)))
        plot(rep(0,length=length(numero)),numero,type="n",axes=F,ylab="",xlab="")
        text(rep(0,length=length(numero)),numero,label=AlphabetAA[numero],col=VectcolAA)
        dev.off()
    }
}



#-----------------------------------------
# Step 3 : compute the number of LS 
# by position
#-----------------------------------------
taille = length(unlist(strsplit(as.character(LS[1,]), "")))
matLS = matrix(NA,ncol = taille, nrow=dim(LS)[1])

for (i in 1:dim(LS)[1]){
   ligne = as.character(LS[i,])
   vect = unlist(strsplit(ligne,""))
   matLS[i,] = vect
}

for (i in 1:dim(matLS)[1]){
   indSel = which(matLS[i,]=="-")
   matLS[i, indSel] = "NA"
}


#recherche le nombre de LS par position
nbrLSbyPos = apply(matLS, 2, f1)
#recherche le nombre de LS dans la même SS par position
nbrLSbyPosbySS = apply(matLS, 2, f2)

#-----------------------------------------
# Step 4 : LS maxtrix
# evolution of the number of LS by position
#-----------------------------------------

Mat3 = matrix(NA, ncol = ncol(matLS), nrow= nrow(matLS))

for (Ibs in 2:length(Alphabet)){
    bs=Alphabet[Ibs]
    Mat3[which(matLS==bs)]=as.numeric(Ibs)
}
rownames(Mat3) = Listpdb


#-----------------------------------------
# Step 5 : LS maxtrix
# visualization of the AA-alignment
#-----------------------------------------

pdf(file.SLalign,paper="USr", width=12)
par(mar = c(3, 3.5,3 ,0))   ###c(bottom, left, top, right). Défaut : 5.1 4.1 4.1 2.1
nf = layout(matrix(c(1,2),1,2,byrow=TRUE), c(22.7,4), c(20,20), TRUE)
image(t(Mat3)[,nrow(Mat3):1],ylab="",cex.lab=1.5,cex.main=1.5,col=VectcolSL,axes=F,main="Structural letter alignment")
      #breaks=seq(min(Mat3,na.rm=T),(max(Mat3,na.rm=T)+1))
axis(2,seq(1,0,length=dim(Mat3)[1]),rownames(Mat3),cex.axis=0.7,las=1)
axis(1,seq(0,(dim(Mat3)[2]-1), by=10)/dim(Mat3)[2],seq(1,dim(Mat3)[2], by=10),cex.axis=0.6,las=1)
#axis(4,seq(1,dim(Mat3)[1], by=5)/dim(Mat3)[1],sort(seq(1,dim(Mat3)[1], by=5),decreasing=T),cex.axis=0.6,las=1)

box()
numero=as.numeric(names(table(Mat3)))
plot(rep(0,length=length(numero)),numero,type="n",axes=F,ylab="",xlab="")
text(rep(0,length=length(numero)),numero,label=Alphabet[numero],col=VectcolSL)	

dev.off()





if (nbgraph>1){
    tailleLim = round(seq(0,dim(Mat3)[1], length=(nbgraph+1)),0)
    for(i in 1:(length(tailleLim)-1)){
        file.SLaligni = paste(pathFig, paste("SL-alignment_", i, ".pdf", sep=""), sep="/")
        pdf(file.SLaligni,width=12,paper="USr")
        par(mar = c(3, 3.5,3 ,0))   ###c(bottom, left, top, right). Défaut : 5.1 4.1 4.1 2.1
        nf = layout(matrix(c(1,2),1,2,byrow=TRUE), c(22.7,4), c(20,20), TRUE)
        ssMat = Mat3[(tailleLim[i]+1):tailleLim[(i+1)], ]

        image(t(ssMat)[,nrow(ssMat):1],ylab="",cex.lab=1.5,cex.main=1.5,col=VectcolSL,axes=F,main="Structural letter alignment")

        axis(2,seq(1,0,length=dim(ssMat)[1]),rownames(ssMat),cex.axis=0.7,las=1)
        axis(1,seq(0,(dim(ssMat)[2]-1), by=10)/dim(ssMat)[2],seq(1,dim(ssMat)[2], by=10),cex.axis=0.6,las=1)

        box()
        numero=as.numeric(names(table(Mat3)))
        plot(rep(0,length=length(numero)),numero,type="n",axes=F,ylab="",xlab="")
        text(rep(0,length=length(numero)),numero,label=Alphabet[numero],col=VectcolSL)	
        dev.off()
       
    }
}




#-----------------------------------------
# Step 6 : Calcul le Neq
#-----------------------------------------
neqLSVect = computNeq(matLS)
neqAAVect = computNeq(matAA)


#-----------------------------------------
# Step 7 : Graphic Representation
# evolution of the number of LS by position
#-----------------------------------------

##def graphics parameters
mm = rbind( neqAAVect, neqLSVect)
valmaxAA = ceiling(max(neqAAVect))
valmaxLS = ceiling(max(neqLSVect))

# coloration according to SS
ssCol2 = f.defColSS(mm) 

#position vector
vcolLeg = rep("white", length=dim(mm)[2])
vcolLeg[seq(1,dim(mm)[2],by=10)]="black"

vLeg = rep("", length=dim(mm)[2])
vLeg[seq(1,dim(mm)[2],by=10)]=seq(1,dim(mm)[2],by=10)




#graphic
pdf(file.neq, paper="USr", width=12)

plot(1:dim(mm)[2],rep(0, length=length(mm[2,])), las=2 ,type="n",ylim=c((-valmaxAA-2),(valmaxLS+2)), axes=FALSE,
     xlab="Alignment positions", las=2, ylab="")

#horizontal and vertical dashed lines
abline(h=c( ((-valmaxAA-1):-1),0,c(1:valmaxLS)), lty=3, col = "gray70")
abline(v= seq(1,dim(mm)[2],by=10), lty=3, col = "gray70")


#plot le Neq LS coloré suivant SS changes
segments(1:dim(mm)[2], rep(1, length=dim(mm)[2]), 1:dim(mm)[2], (mm[2,]+1),col=ssCol2,lwd=2, lty=1)
#plot le Neq AA 
segments(1:dim(mm)[2],rep(-1, length=dim(mm)[2]), 1:dim(mm)[2], (-mm[1,]-1),lwd=2, col = "orange",lty=1)
axis(2,((-valmaxAA-2):(valmaxLS+2)), c(NA,abs(-valmaxAA:0),NA,(0:(valmaxLS)),NA),las=2)
#legend of NeqSL
legend("topleft", horiz = TRUE, legend=c("helix", "strand","loop","SS changes"), 
       col = c("red","magenta","gray","purple"), cex=0.8, pch=15,x.intersp=1, title=expression(neq[SL]))
#position vector
points(1:dim(mm)[2],rep(0, length=dim(mm)[2]), pch=124, col = vcolLeg)
text(1:dim(mm)[2],rep(-0.9, length=dim(mm)[2]),vLeg , col = 1, srt=90, cex=0.5, pos=3, offset=0.8)
#ylab of Neq
mtext(expression(neq[SL]),side=2,line=2, srt=90, at=2.5)
mtext(expression(neq[AA]),side=2,line=2, srt=90, at=-1.7,col="orange")

dev.off()


#-----------------------------------------
# Step 8 : file output of all positions
#-----------------------------------------
ssVect = apply(matLS,2,f.varSS)
nbrGapVectAA =  f.nbrGap(matAA)
nbrGapVectLS =  f.nbrGap(matLS)


signifVect = rep(" ", length=length(neqAAVect))
signifVect[which(neqAAVect>1)] = "!"
signifVect[which(neqLSVect>seuilVar)] = "/"
signifVect[which( (neqLSVect>seuilVar) & (nbrLSbyPosbySS> 1) )] = "*"
tmp = which( (neqLSVect>seuilVar) & (nbrLSbyPosbySS> 1) )
tmp2 = which(ssVect=="v")
signifVect[intersect(tmp,tmp2)] = "**"





matoutput = cbind(data.frame(fileC , f.ExtractPDBRef.Pos(1:length(nbrAAbyPos), folder_out), 
                             nbrAAbyPos, nbrLSbyPos,round(neqAAVect,2),
                                round(neqLSVect,2)), cbind(ssVect,signifVect ))
colnames(matoutput) = c("aligned.pos",colnames(fileC)[-1],paste("pdb.",f.ExtractPDBRef(folder_out),".pos",sep=""),"nbr.AA","nbr.SL","neqAA","neqSL","SS","signif")


write.table(matoutput, file=file.out,quote=F, col.names=T, row.names=F, sep=";")


#-----------------------------------------
# Step 9: Determination of the important positions
#-----------------------------------------

ph0.0 = paste("********* The MSA has", length(neqAAVect>=1), "aligned positions *********\n",sep=" ")
write(ph0.0, file = file_out1)



write("**************************************************\n*** Occurrences of different type of positions ***\n**************************************************\n",
     file = file_out1, append=T)


#1- position conservées en termes de AA
m1.AA = round(mean(neqAAVect[which(neqAAVect>0)]),2)
s1.AA = round(sd(neqAAVect[which(neqAAVect>0)]),2)

ph0 = paste("* ", length(which(neqAAVect>=1)), " positions with a computed neq_AA value: average neq_AA = ", 
            m1.AA, " (+/-", s1.AA, ")",sep="")
write(ph0, file = file_out1, append=T)

posAAConserv = which(neqAAVect==1)
posAAmut = which(neqAAVect>1)
posAAMoyVar =intersect(which(neqAAVect>1), which(neqAAVect<seuilVar))
posAAVar = which(neqAAVect>=seuilVar)

pos.type.AA = c(length(posAAConserv),length(posAAmut))

ph1 = paste("    -AA conserved positions : ", length(posAAConserv), " (= ",  round(length(posAAConserv)/sum(pos.type.AA),3)*100 ,"% of positions)", sep=" ")
ph2.1 = paste("    -AA weakly variable positions : ", length(posAAMoyVar), " (= ",  round(length(posAAMoyVar)/sum(pos.type.AA),3)*100 ,"% of positions)", sep=" ")
ph2.2 = paste("    -AA strongly variable positions : ", length(posAAVar), " (= ",  round(length(posAAVar)/sum(pos.type.AA),3)*100 ,"% of positions)", sep=" ")


write(paste(ph1,ph2.1,ph2.2,sep="\n"), file = file_out1, append=T)





#2- position conservées en termes de SL
m1.SL = round(mean(neqLSVect[which(neqLSVect>0)]),2)
s1.SL = round(sd(neqLSVect[which(neqLSVect>0)]),2)

ph01 = paste("\n* ", length(which(neqLSVect>=1)), " positions with a computed neq_SL value: average neq_SL = ", 
            m1.SL, " (+/-", s1.SL, ")",sep="")
write(ph01, file = file_out1, append=T)


posSLConserv = which(neqLSVect==1)
posSLMoyVar = intersect(which(neqLSVect>1), which(neqLSVect<seuilVar))
posSLVar = which(neqLSVect>=seuilVar)
pos.type.SL = c(length(posSLConserv) , length(posSLMoyVar) , length(posSLVar))

ph3 = paste("    -SL conserved positions : ", length(posSLConserv), " (= ",  round(length(posSLConserv)/sum(pos.type.SL),3)*100 ,"% of positions)", sep=" ")
ph4 = paste("    -SL weakly variable positions : ", length(posSLMoyVar), " (= ",  round(length(posSLMoyVar)/sum(pos.type.SL),3)*100 ,"% of positions)", sep=" ")
ph5 = paste("    -SL strongly variable positions : ", length(posSLVar), " (= ",  round(length(posSLVar)/sum(pos.type.SL),3)*100 ,"% of positions)", sep=" ")

write(paste(ph3,ph4,ph5,sep="\n"), file = file_out1, append=T)


#-----------------------------------------
# Step 10 : Analysis of the repartition of the different positions along the sequence
#-----------------------------------------
write("\n\n****************************************************\n*** Distribution of position types along the MSA ***\n****************************************************",
      file = file_out1, append=T)

f.WriteConservAlongSeq(neqLSVect, type="Conser", file = file_out1)
f.WriteConservAlongSeq(neqLSVect, type="WeakVar", file = file_out1)
f.WriteConservAlongSeq(neqLSVect, type="StrongVar", file = file_out1)

#-----------------------------------------
# Step 11: mutation analysis
#-----------------------------------------

if (length(posAAmut) != 0){
  for (posi in 1:length(posAAmut)){
    pos = posAAmut[posi]
	
    ph0 = paste("Mutation n°",posi, "- aligned position:", pos,sep=" ")
    ph1 = paste("*** Repartition of AA in this mutated position ***")

    write(paste(ph0,ph1,sep="\n"), file = file.mut,append=T)
    tc = table(matAA[,pos])
    if(is.element("NA", names(tc))=="TRUE") { 
    	idx = which(names(tc)=="NA")
    	names(tc)[idx] = "-"
    }
	for(i in 1:length(tc)){
		ph2 = paste("    ", names(tc)[i], ": ", tc[i],sep=" ")
		write(ph2, file = file.mut, append=T)
	}
	ph3 = paste("*** Repartition of SL in this mutated position ***")  
	write(ph3, file = file.mut, append=T)
	for(j in 1:length(tc)){
		aa = names(tc)[j]
		
		if (aa != "-"){
		    posSel = which(matAA[,pos]==aa)
		    tc2 = sort(table(matLS[posSel,pos]),decreasing=T)
		    if(is.element("NA", names(tc2))=="TRUE") { 
       	               idx = which(names(tc2)=="NA")
    	                names(tc2)[idx] = "-"
                    }
		    
		    p0 = paste("    ",aa, "(", tc[j],")     ", sep="")
                    p1 = paste(paste(names(tc2),tc2,sep=":"), collapse=" ; ")
                    write(paste(p0, p1, sep=" "), file = file.mut, append=T)
		}else{
		    p0 = paste("    ",aa, "(", tc[j],")     ", sep="")
            p1 = 0
            write(paste(p0, p1, sep=" "), file = file.mut, append=T)
        }
	}		
	write("--------------------------------------------------------", file = file.mut, append=T)
  }
}

if (length(posAAmut) == 0){
   print ("zero mutation occurs in the sequences")
}





#-----------------------------------------
# Step 11'': Structural variable positions analysis
#-----------------------------------------
posSLVar12 = sort(c(posSLVar,posSLMoyVar))

if (length(posSLVar12) != 0){
  for (posi in 1:length(posSLVar12)){
    pos = posSLVar12[posi]
	
    ph0 = paste("Structural variable position n°", posi, "- aligned position:", pos,"neqSL=",round(neqLSVect[pos],2),sep=" ")
    ph1 = paste("*** Repartition of SL in this structural variable position ***")

    write(paste(ph0,ph1,sep="\n"), file = file.SLVar,append=T)
    tc = table(matLS[,pos])
    if(is.element("NA", names(tc))=="TRUE") { 
    	idx = which(names(tc)=="NA")
    	names(tc)[idx] = "-"
    }
	for(i in 1:length(tc)){
		ph2 = paste("    ", names(tc)[i], ": ", tc[i],sep=" ")
		write(ph2, file = file.SLVar, append=T)
	}
		
	write("--------------------------------------------------------------", file = file.SLVar, append=T)
  }
}

if (length(posSLVar12) == 0){
   print ("there is no structural variable positions in the MSLA")
}




#-----------------------------------------
# Step 12: Create file with the 
# correspondence between 
# - aligned position
# - uniProt position
# - position in pdb in all pdb files
#-----------------------------------------
list.pos = f.ExtractPDB.Pos(folder_out)

mat.all.pos = cbind(fileC, list.pos[[1]])
colnames(mat.all.pos) = c(colnames(fileC), list.pos[[2]])

write.table(mat.all.pos, file=filePos_out,quote=T, col.names=T, row.names=F, sep=";")


