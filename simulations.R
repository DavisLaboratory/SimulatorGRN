#----sim 1----
net.edge=matrix(c('A','act','C','B','rep','D','C,D','act,rep','E','E','act','C','A,B','rep,rep','F'),byrow=T,ncol=3)
net.edge=as.data.frame(net.edge,stringsAsFactors=F)
colnames(net.edge)=c('regulator','interaction','target')
net.nodes=data.frame('gene'=c('A','B','C','D','E','F'),stringsAsFactors=F)
grn1=list('nodes'=net.nodes,'edges'=net.edge)

grn=model.init(grn1)
igrn=grn.to.igrn(grn)
ext=c(0.5,0.3)
names(ext)=LETTERS[1:2]
res=model.solve(igrn,ext)

#----sim 2----
net.edge=matrix(c('A','act','D','',
				  'A','act','L','',
				  'B','act','E','',
				  'C','act','F','',
				  'D','act','G','',
				  'E','act','I','',
				  'F','act','K','',
				  'K','act','J','',
				  'D,E','act,act','H','and',
				  'E,F','act,rep','J','and')
				,byrow=T,ncol=4)
net.edge=as.data.frame(net.edge,stringsAsFactors=F)
colnames(net.edge)=c('regulator','interaction','target','type')
net.nodes=data.frame('gene'=LETTERS[1:12],stringsAsFactors=F)
grn2=list('nodes'=net.nodes,'edges'=net.edge)

grn=model.init(grn2)
igrn=grn.to.igrn(grn)
ext=c(0.5,0.3,0.7)
names(ext)=LETTERS[1:3]
res=model.solve(igrn,ext)
