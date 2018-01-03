

require(igraph)
require(bipartite)
source('utils.r')
#### Step 1: Read input data.
#### In this example case the interactions are reported as adjacency lists per area

input_data <- read.csv('./DATA_NAR/Canada/Bottomup_Lakeweb.csv', stringsAsFactors = F)

#### Step 2: From this table we create networks for each one of the local communities

areas <- sort(unique(input_data$Lake_ID))

#### Step 3: If data has to be aggregated set this variable to TRUE
aggregate <- FALSE
if(aggregate){
  whole_net <- empty_graph()
}

bipartite <- FALSE
quantified <- FALSE

species <- c()
links <- c()
connectances <- c()
links_per_sp <- c()
indegree <- c()
indegree_normalized <- c()
outdegree <- c()
outdegree_normalized <- c()
sd_vul <- c()
sd_gen <- c()
basal <- c()
top <- c()
intermediate <- c()
S_basal <- c()
S_intermediate <- c()
S_top <- c()
omnivory <- c()
mfcls <- c()

#### this is the output table where the results are going to be stored
output <- NULL
for(a in areas){
  cur_com <- input_data[which(input_data$Lake_ID == a),]
  
  cur_com[is.na(cur_com)] <- 0
  row.names(cur_com) <- cur_com$Prey
  n <- as.matrix(cur_com[-c(1:3)])
  
  remove <- which(rowSums(n)==0 & colSums(n)==0)
  if(length(remove) != 0){
    n <- n[-remove, -remove]
  }
  local_net <- graph_from_adjacency_matrix(n)
  
  
  # local_net <- make_empty_graph()
  # for(int in 1:dim(cur_com)[1]){
  #   prey <- cur_com[int,]$prey
  #   if(! prey %in% V(local_net)$name){
  #     local_net <- local_net + vertices(prey)
  #   }
  #   pred <- cur_com[int,]$pred
  #   if(! pred %in% V(local_net)$name){
  #     local_net <- local_net + vertices(pred)
  #   }
  #   local_net[prey, pred] <- 1
  #   if(quantified){
  #     E(local_net)[prey,pred]$strength <- cur_com[int,]$strength
  #   }
  # }
  
  if(aggregate){
    local_net <- union(whole_net, local_net)
    whole_net <- local_net
  }

  #### Step 4: Calculate the network features for the current community
  if(!bipartite){
    S <- vcount(local_net)
    species <- append(species, S)
    ls <- ecount(local_net)
    links <- append(links, ls)
    connectances <- append(connectances, (2*ls/((S-1)*S)))
    links_per_sp <- append(links_per_sp, ls/S)
    
    n <- as_adjacency_matrix(local_net)
   
    indegree <- append(indegree, MeanGenerality(n))
    indegree_normalized <- append(indegree_normalized, NormalisedGenerality(n))
    outdegree <- append(outdegree, MeanVulnerability(n))
    outdegree_normalized <- append(outdegree_normalized, NormalisedVulnerability(n))
    
    sd_gen <- append(sd_gen, SDGenerality(n))
    sd_vul <- append(sd_vul, SDVulnerability(n))
    
    #####
    basal <- append(basal, FractionOfBasal(n))
    top <- append(top, FractionOfTop(n))
    intermediate <- append(intermediate, FractionOfIntermediate(n))
    S_basal <- append(S_basal, NumberOfBasal(n))
    S_top <- append(S_top, NumberOfTop(n))
    S_intermediate <- append(S_intermediate, NumberOfIntermediate(n))
    
  
    omnivory <- append(omnivory, Omnivory(n))
    
    mfcl <- tryCatch({
      MeanFoodChainLength(n)
    }, warning = function(w) {
      'NA'
    }, error = function(e) {
      'NA'
    }, finally = {

    })
    mfcls <- append(mfcls, mfcl)
  }else{
    #### if the network is bipartite we convert it to the bipartite representation
    g_bipart <- graph.bipartite(bipartite.mapping(local_net)$type, as.vector(t(get.edges(local_net, 1:length(E(local_net))))), directed=T)
    V(g_bipart)$name <- V(local_net)$name
    
    S <- vcount(local_net)
    Sc <- length(which(V(gbipart)$type == TRUE))
    Sr <- length(which(V(gbipart)$type == FALSE))
    
    n <- as_incidence_matrix(g_bipart)
    
    species <- append(species, S)
    ls <- ecount(local_net)
    links <- append(links, ls)
    connectances <- append(connectances, networklevel(n, 'connectance'))
    links_per_sp <- append(links_per_sp, networklevel(n, 'links per species'))
    
    indegree <- append(indegree, networklevel(n, 'generality'))
    #indegree_normalized <- append(indegree_normalized, NormalisedGenerality(n))
    outdegree <- append(outdegree, networklevel(n, 'vulnerability'))
    #outdegree_normalized <- append(outdegree_normalized, NormalisedVulnerability(n))
    
    sd_gen <- append(sd_gen, SDGenerality(n))
    sd_vul <- append(sd_vul, SDVulnerability(n))
  }
  
  # cur_out <- data.frame(areas, species, links, connectances, links_per_sp, indegree) #, indegree_normalized, outdegree, outdegree_normalized, sd_gen, sd_vul, basal, top, intermediate, S_basal, S_top, S_intermediate, omnivory, mfcls)
  # 
  # if(is.null(output)){
  #   output <- cur_out
  # }else{
  #   output <- rbind(output, cur_out)
  # } 
}


metadata <- read.csv('./DATA_NAR/Canada/Lake_Environment128.csv', stringsAsFactors = F)
areas_ha <- metadata[match(areas, metadata$Wby_LID),]$Area_ha
cur_out_bottom <- data.frame(areas=areas_ha, species, links, connectances, links_per_sp, indegree, outdegree, sd_gen, sd_vul, basal, top, intermediate, S_basal, S_top, S_intermediate, omnivory, mfcls) 

cur_out$type <- 'TD'
cur_out_bottom$type <- 'BU'

cur_out <- rbind(cur_out, cur_out_bottom)

cur_out <- cur_out[-which(is.na(cur_out$areas)),]

require(ggplot2)

ggplot(cur_out, aes(log10(areas), outdegree, colour=type)) + theme_bw() + geom_jitter() + stat_smooth(method = "lm")

