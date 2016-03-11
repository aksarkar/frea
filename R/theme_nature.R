requireNamespace('ggplot2')
requireNamespace('grid')

#' @importFrom grid unit

theme_nature <- theme(
    line=element_line(color='black', size=(0.5 / ggplot2:::.pt), linetype=1,
                      lineend='square'),
    rect=element_blank(),
    text=element_text(family='Helvetica', face='plain', color='black', size=7,
                      hjust=0.5, vjust=0.5, angle=0, lineheight=1,
                      margin=unit(rep(2, 4), 'points'), debug=FALSE),
    axis.line=element_line(),
    axis.line.x=element_line(),
    axis.line.y=element_line(),
    axis.text=element_text(size=5),
    axis.ticks=element_line(color='black'),
    axis.ticks.length=unit(2, 'points'),
    axis.title=element_text(margin=unit(rep(4, 4), 'points')),
    axis.title.x=element_text(),
    axis.title.y=element_text(angle=90),
    legend.background=element_blank(),
    legend.box=NULL,
    legend.direction=NULL,
    legend.justification='center', 
    legend.key.height=NULL,
    legend.key.size=unit(7, 'points'),
    legend.key.width=NULL,
    legend.key=element_rect(fill='grey95', color='white'),
    legend.margin=unit(0.2, 'cm'),
    legend.position='none',
    legend.text.align=NULL,
    legend.text=element_text(size=rel(0.8)),
    legend.title.align=NULL, 
    legend.title=element_text(size=rel(0.8), face='bold', hjust=0),
    panel.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank(),
    panel.margin=unit(0, 'points'),
    plot.background=element_blank(),
    plot.margin=unit(rep(0, 4), 'points'),
    plot.title=element_blank(),
    strip.background=element_blank(),
    strip.text.x=element_text(),
    strip.text.y=element_text(angle=-90),
    strip.text=element_text(),
    complete=TRUE)

line_attrs <- list(size=I(.35 / ggplot2:::.pt))

update_geom_defaults('line', line_attrs)
update_geom_defaults('path', line_attrs)
update_geom_defaults('abline', line_attrs)
update_geom_defaults('point', list(size=.5))
