# This script basically plots networks using an pre-build igraph object from RTN analysis. 
# The functions only work for AMAP and RMAP. AMAPDEND still not implemented.

# graphAmap() function ####
graphAmap <- function(g,dat,refcol=1,alias,pvalue)
{
  rdp <- RedeR::RedPort()
  gr <- RedeR::att.mapv(g=g, dat=dat, refcol=1)
  gr <- RedeR::att.setv(g=gr, from=colnames(dat)[alias], to="nodeAlias")
  gr <- RedeR::att.setv(g=gr, from=colnames(dat)[pvalue], to="nodeColor", 
                        breaks=seq(min(dat[,pvalue]), max(dat[,pvalue]),0.4), pal=2)
  RedeR::calld(rdp)
  RedeR::addGraph(rdp, gr)
  RedeR::relax(rdp)
  scl <- gr$legNodeColor$scale
  leg <- gr$legNodeColor$legend
  RedeR::addLegend.color(rdp, colvec=scl, labvec=leg, title="node color (p-value)")
  scl <- gr$legNodeSize$scale
  leg <- as.character(round(as.numeric(gr$legNodeSize$legend)))
  RedeR::addLegend.size(rdp, sizevec=scl, labvec=leg, title="node size (Regulon Size)",
                        position="bottomright", intersp=10,ftsize=10,vertical=F)
  scl <- gr$legEdgeWidth$scale
  leg <- gr$legEdgeWidth$legend
  RedeR::addLegend.size(rdp,type="edge", sizevec=scl, labvec=leg, title="edge width (Weigth)",
                        position="topleft", intersp=10,ftsize=10,vertical=T)
}

# graphAmap arguments:
# g = graph object (amap)
# dat = pheno data.frame
# refcol = column within pheno data with probe ids
# alias = string with name of respective column (node symbol)
# pvalue = same for alias, but with p-values

# graphRmap() function ####
graphRmap <- function(g,dat,refcol=1,colid=3)
{
  rdp <- RedeR::RedPort()
  gr <- RedeR::att.mapv(g=g, dat=dat, refcol=1)
  gr <- RedeR::att.setv(g=gr, from=colnames(dat)[colid], to="nodeAlias")
  RedeR::calld(rdp)
  RedeR::addGraph(rdp, gr)
  RedeR::relax(rdp)
  scl <- gr$legEdgeColor$scale
  leg <- gr$legEdgeColor$legend
  RedeR::addLegend.color(rdp,type="edge", colvec=scl,labvec=leg,title="TF regulation", position="topright",vertical=T)
}

# graphRmap arguments:
# g = graph object (rmap)
# dat = pheno data.frame
# refcol = column within pheno data with probe ids
# Column 3 MUST hold SYMBOL. Otherwise,set other column. Default 3.