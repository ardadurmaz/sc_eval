library(Matrix)
library(ggplot2)

suppressMessages(source('~/sc_eval/scripts/sc_utils.R'))
suppressMessages(source('~/sc_eval/scripts/sc_functions.R'))
set.seed(42)

#uniq.labels <- c('H3122-1', 'H3122-2', 'H3122-Alec-4h', 'H3122-Alec-48h', 'H3122-Alec-3W', 'H3122-erAlec-1', 'H3122-GFP-erAlec',
                 #'H3122-Lor-48h', 'H3122-Lor-3W', 'H3122-erLor-1', 'H3122-erLor',
                 #'H3122-Cz-48h', 'H3122-Cz-3W', 'H3122-erCz-3', 'H3122-erCz')
#uniq.colors <- setNames(ggsci::pal_simpsons()(15), uniq.labels)

#uniq.labels <- c('0h', '12h', '3h', '6h')
#uniq.colors <- setNames(ggsci::pal_jco()(4), uniq.labels)
uniq.labels <- c('Week-0', 'Week-1', 'Week-2', 'Week-6')
uniq.colors <- setNames(ggsci::pal_jco()(4), uniq.labels)

title_sub <- 'EGEOD103334'

for(i in 1:7){
  ## Slingshot ##
  title <- paste0(title_sub, sprintf('_S%d', i))
  setwd(sprintf('~/sc_eval/data/%s', title))
  for(d_type in c('raw', 'scimpute')){
    for(n_type in c('deconv', 'sctransform', 'dca')){
      for(r_type in c('umap', 'tsne', 'dhaka', 'dm', 'paga_umap')){
        if(d_type == 'scimpute' && n_type == 'dca')
          next
        
        tag <- sprintf('%s_%s_%s_%s_Slingshot.rds', title, d_type, n_type, r_type)
        traj_res <- readRDS(tag)[[1]]
        #traj_res$Label <- factor(traj_res$Label, levels = c('H3122-1', 'H3122-2', 'H3122-Alec-4h', 'H3122-Alec-48h', 'H3122-Alec-3W', 'H3122-erAlec-1', 'H3122-GFP-erAlec'))
        #traj_res$Label <- factor(traj_res$Label, levels = c('H3122-1', 'H3122-2', 'H3122-Lor-48h', 'H3122-Lor-3W', 'H3122-erLor', 'H3122-erLor-1'))
        #traj_res$Label <- factor(traj_res$Label, levels = c('H3122-1', 'H3122-2', 'H3122-Cz-48h', 'H3122-Cz-3W', 'H3122-erCz', 'H3122-erCz-3'))
        #traj_res$Label <- factor(traj_res$Label, levels = c('0h', '3h', '6h', '12h'))
        traj_res$Label <- factor(traj_res$Label, levels = c('Week-0', 'Week-1', 'Week-2', 'Week-6'))
        
        p <- ggplot() +
          geom_point(subset(traj_res, traj_res$Type=='CB'), mapping=aes(x=Dim.1, y=Dim.2, color=Label), shape=20, size=1) +
          theme_classic() +
          labs(color='') +
          xlab('Dimension 1') +
          ylab('Dimension 2') +
          theme(axis.title = element_text(face='bold', size=12)) +
          #scale_color_manual(values=uniq.colors, labels=setNames(c('Naive-1', 'Naive-2', 'Alec., 4h', 'Alec., 48h', 'Alec., 3w', 'erAlec., 1', 'erAlec., 2'), uniq.labels[1:7])) +
          #scale_color_manual(values=uniq.colors, labels=setNames(c('Naive-1', 'Naive-2', 'Lor., 48h', 'Lor., 3w', 'erLor., 1', 'erLor., 2'), uniq.labels[c(1,2,8,9,10,11)])) +
          #scale_color_manual(values=uniq.colors, labels=setNames(c('Naive-1', 'Naive-2', 'Criz., 48h', 'Criz., 3w', 'erCriz., 1', 'erCriz., 2'), uniq.labels[c(1,2,12,13,14,15)])) +
          scale_color_manual(values=uniq.colors) +
          guides(color=guide_legend(override.aes = list(size=4)))
        
        for(s in setdiff(unique(traj_res$Type), 'CB')){
          p <- p + geom_path(subset(traj_res, traj_res$Type==s), mapping=aes(x=Dim.1,y=Dim.2), color='black')
        }
        
        ggsave(p, filename = sprintf('%s_%s_%s_%s_Slingshot.pdf', title, d_type, n_type, r_type), width = 6.5, height = 4)
      }
    }
  }
}

for(i in 1:7){
  ## Monocle ##
  title <- paste0(title_sub, sprintf('_S%d', i))
  setwd(sprintf('~/sc_eval/data/%s', title))
  for(d_type in c('raw', 'scimpute')){
    for(n_type in c('deconv', 'sctransform', 'dca')){
      if(d_type == 'scimpute' && n_type == 'dca')
        next
      
      tag <- sprintf('%s_%s_%s_Monocle2.rds', title, d_type, n_type)
      traj_res <- readRDS(tag)[[1]]
      #traj_res$Label <- factor(traj_res$Label, levels = c('H3122-1', 'H3122-2', 'H3122-Alec-4h', 'H3122-Alec-48h', 'H3122-Alec-3W', 'H3122-erAlec-1', 'H3122-GFP-erAlec'))
      #traj_res$Label <- factor(traj_res$Label, levels = c('H3122-1', 'H3122-2', 'H3122-Lor-48h', 'H3122-Lor-3W', 'H3122-erLor', 'H3122-erLor-1'))
      #traj_res$Label <- factor(traj_res$Label, levels = c('H3122-1', 'H3122-2', 'H3122-Cz-48h', 'H3122-Cz-3W', 'H3122-erCz', 'H3122-erCz-3'))
      #traj_res$Label <- factor(traj_res$Label, levels = c('0h', '3h', '6h', '12h'))
      traj_res$Label <- factor(traj_res$Label, levels = c('Week-0', 'Week-1', 'Week-2', 'Week-6'))
      
      p <- ggplot() +
        geom_point(subset(traj_res, !is.na(traj_res$Label)), mapping=aes(x=Dim.1, y=Dim.2, color=Label), shape=20, size=1.5, stroke=0) +
        geom_point(subset(traj_res, is.na(traj_res$Label)), mapping=aes(x=Dim.1, y=Dim.2), color='black', shape=20, size=0.75, stroke=0, alpha=0.8) +
        theme_classic() +
        labs(color='') +
        xlab('Dimension 1') +
        ylab('Dimension 2') +
        theme(axis.title = element_text(face='bold', size=12)) +
        #scale_color_manual(values=uniq.colors, labels=setNames(c('Naive-1', 'Naive-2', 'Alec., 4h', 'Alec., 48h', 'Alec., 3w', 'erAlec., 1', 'erAlec., 2'), uniq.labels[1:7])) +
        #scale_color_manual(values=uniq.colors, labels=setNames(c('Naive-1', 'Naive-2', 'Lor., 48h', 'Lor., 3w', 'erLor., 1', 'erLor., 2'), uniq.labels[c(1,2,8,9,10,11)])) +
        #scale_color_manual(values=uniq.colors, labels=setNames(c('Naive-1', 'Naive-2', 'Criz., 48h', 'Criz., 3w', 'erCriz., 1', 'erCriz., 2'), uniq.labels[c(1,2,12,13,14,15)])) +
        scale_color_manual(values=uniq.colors) +
        guides(color=guide_legend(override.aes = list(size=4)))
      ggsave(p, filename = sprintf('%s_%s_%s_Monocle2.pdf', title, d_type, n_type), width = 6.5, height = 4)
    }
  }
}
