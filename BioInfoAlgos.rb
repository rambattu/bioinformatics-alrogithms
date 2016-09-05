class BioInfoAlgos 

    def rna_to_codon_hash
        {
            :AAA => "K", :AAC => "N", :AAG => "K", :AAU => "N", :ACA => "T", :ACC => "T", :ACG => "T", :ACU => "T",
            :AGA => "R", :AGC => "S", :AGG => "R", :AGU => "S", :AUA => "I", :AUC => "I", :AUG => "M", :AUU => "I",
            :CAC => "H", :CAG => "Q", :CAU => "H", :CCA => "P", :CCC => "P", :CCG => "P", :CCU => "P", :CGA => "R",
            :CGC => "R", :CGG => "R", :CGU => "R", :CUA => "L", :CUC => "L", :CUG => "L", :CUU => "L", :GAA => "E",
            :GAC => "D", :GAG => "E", :GAU => "D", :GCA => "A", :GCC => "A", :GCG => "A", :GCU => "A", :GGA => "G",
            :GGC => "G", :GGG => "G", :GGU => "G", :GUA => "V", :GUC => "V", :GUG => "V", :GUU => "V", :UAA => "",
            :UAC => "Y", :UAG => "",  :UAU => "Y", :UCA => "S", :UCC => "S", :UCG => "S", :UCU => "S", :UGA => "",
            :UGC => "C", :UGG => "W", :UGU => "C", :UUA => "L", :UUC => "F", :UUG => "L", :UUU => "F", :CAA => "Q",
        }
    end # end of rna_to_codon_hash function

    def dna_to_rna(dna)
        dna.gsub("T", "U")
    end # end of dna_to_rna function

    def rna_to_dna(rna)
        rna.gsub("U", "T")
    end # end of rna_to_dna function

    def amino_acid_mass_hash
        {
            "G" => 57, "A" => 71, "S" => 87, "P" => 97, "V" => 99,
            "T" => 101, "C" => 103, "I" => 113, "L" => 113, "N" => 114,
            "D" => 115, "K" => 128, "Q" => 128, "E" => 129, "M" => 131,
            "H" => 137, "F" => 147, "R" => 156, "Y" => 163, "W" => 186,
        }
    end # end of amino_acid_mass_hash

    def kmers(text,k)
        # Given a text, give out all the possible kmers in that text
        kmer_a = []
        (0..(text.length-k)).each {|i| kmer_a << text.slice(i,k)}
        return kmer_a
    end # end of kmers function

    def kmer_prefix(kmer)
        prefix = kmer.slice(0, kmer.length-1)
        return prefix
    end

    def kmer_suffix(kmer)
        suffix = kmer.slice(1, kmer.length-1)
        return suffix
    end

    def paired_kmer_prefix(kmer, k, d)
        kmers = kmer.split("|")
        return (kmer_prefix(kmers[0]) + "|" + kmer_prefix(kmers[1]))
    end # end of paired_kmer_prefix

    def paired_kmer_suffix(kmer, k, d)
        kmers = kmer.split("|")
        return (kmer_suffix(kmers[0]) + "|" + kmer_suffix(kmers[1]))
    end # end of paired_kmer_suffix

    def blosum_62 
        lines = IO.readlines("BLOSUM62.txt") 
        #    A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y
        # A  4  0 -2 -1 -2  0 -2 -1 -1 -1 -1 -2 -1 -1 -1  1  0  0 -3 -2
        # C  0  9 -3 -4 -2 -3 -3 -1 -3 -1 -1 -3 -3 -3 -3 -1 -1 -1 -2 -2
        # D -2 -3  6  2 -3 -1 -1 -3 -1 -4 -3  1 -1  0 -2  0 -1 -3 -4 -3
        # Lets get rid of the first line
        lines.shift
        # puts lines
        blosum = {}
        protein_a = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
        i = 0
        lines.each do |line|
            score_a = line.split(/\s+/)
            # puts score_a.join(",")
            # lets get rid of the first entry, that will be letter of protein A, C, D,..
            score_a.shift
            blosum[protein_a[i]] = {
                "A" => score_a[0].to_i, "C" => score_a[1].to_i, "D" => score_a[2].to_i, "E" => score_a[3].to_i, "F" => score_a[4].to_i,
                "G" => score_a[5].to_i, "H" => score_a[6].to_i, "I" => score_a[7].to_i, "K" => score_a[8].to_i, "L" => score_a[9].to_i,
                "M" => score_a[10].to_i, "N" => score_a[11].to_i, "P" => score_a[12].to_i, "Q" => score_a[13].to_i, "R" => score_a[14].to_i,
                "S" => score_a[15].to_i, "T" => score_a[16].to_i, "V" => score_a[17].to_i, "W" => score_a[18].to_i, "Y" => score_a[19].to_i,
            }
            i += 1
        end       

        return blosum
    end # end of blosum_62 

    def pam_250 
        lines = IO.readlines("PAM250_1.txt") 
        #    A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y
        # A  4  0 -2 -1 -2  0 -2 -1 -1 -1 -1 -2 -1 -1 -1  1  0  0 -3 -2
        # C  0  9 -3 -4 -2 -3 -3 -1 -3 -1 -1 -3 -3 -3 -3 -1 -1 -1 -2 -2
        # D -2 -3  6  2 -3 -1 -1 -3 -1 -4 -3  1 -1  0 -2  0 -1 -3 -4 -3
        # Lets get rid of the first line
        lines.shift
        # puts lines
        pam = {}
        protein_a = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
        i = 0
        lines.each do |line|
            score_a = line.split(/\s+/)
            # puts score_a.join(",")
            # lets get rid of the first entry, that will be letter of protein A, C, D,..
            score_a.shift
            pam[protein_a[i]] = {
                "A" => score_a[0].to_i, "C" => score_a[1].to_i, "D" => score_a[2].to_i, "E" => score_a[3].to_i, "F" => score_a[4].to_i,
                "G" => score_a[5].to_i, "H" => score_a[6].to_i, "I" => score_a[7].to_i, "K" => score_a[8].to_i, "L" => score_a[9].to_i,
                "M" => score_a[10].to_i, "N" => score_a[11].to_i, "P" => score_a[12].to_i, "Q" => score_a[13].to_i, "R" => score_a[14].to_i,
                "S" => score_a[15].to_i, "T" => score_a[16].to_i, "V" => score_a[17].to_i, "W" => score_a[18].to_i, "Y" => score_a[19].to_i,
            }
            i += 1
        end       

        return pam
    end # end of PAM_250 


    def shared_kmers(k, x, y)
        x_kmers = kmers(x, k)
        y_kmers = kmers(y, k)

        # puts "Obtained kmers"

        pairs = []

        x_kmer_h = {}
        x_kmers.each_index do |idx|
            kmer = x_kmers[idx]
            unless x_kmer_h.has_key?(kmer)
                x_kmer_h[kmer] = [idx]                 
            else
                x_kmer_h[kmer] << idx
            end 

            kmer = reverse_complement(x_kmers[idx])
            unless x_kmer_h.has_key?(kmer)
                x_kmer_h[kmer] = [idx]                 
            else
                x_kmer_h[kmer] << idx
            end 


        end

        # puts "Finished pouplating x_kmer hash"

        y_kmer_h = {}
        y_kmers.each_index do |idx|
            kmer = y_kmers[idx]
            unless y_kmer_h.has_key?(kmer)
                y_kmer_h[kmer] = [idx]                 
            else
                y_kmer_h[kmer] << idx
            end 
        end

        # puts "Finished pouplating y_kmer hash"

        x_kmer_h.keys.each do |x_kmer|
            x_kmer_h[x_kmer].each do |x_idx|
                if y_kmer_h.has_key?(x_kmer)
                    y_kmer_h[x_kmer].each do |y_idx|
                        pairs << [x_idx, y_idx]
                    end
                end
            end
        end


        pairs.sort! do |a, b|
            a[1] <=> b[1]
            # a[0] <=> b[0]
        end
        return pairs
    end # end of shared_kmers function


    def graph_to_genome(colored_edges)
        <<-DOC
        GraphToGenome(GenomeGraph)
             P ← an empty set of chromosomes
             for each cycle Nodes in GenomeGraph
                  Chromosome ← CycleToChromosome(Nodes)
                  add Chromosome to P
             return P


        Sample Input:
        (2, 4), (3, 6), (5, 1), (7, 9), (10, 12), (11, 8)
        Sample Output:
        (+1 -2 -3)(-4 +5 -6)        
        DOC

        # We have to detect cycle first and then form chromosme
        # The cycle ends when a second digit in the edge is less than the first digit

        p = []
        list = []
        lowest = -1
        colored_edges.each do |edge|
            if lowest == -1
                list = edge
                lowest = edge.min 
            elsif (edge.min < lowest)
                list << edge[0]
                list.unshift(edge[1])
                p << list
                list = []
                lowest = -1
            else
                list << edge[0]
                list << edge[1]
            end
        end
        # p << list unless list.empty?
        if list.length >= 1
            last_elem = list.pop  
            list.unshift(last_elem)
            p << list
        end

        genome = []
        p.each do |chrome|
            # puts chrome.join(" ")
            # puts ""
            # puts cycle_to_chromosome(chrome).join(" ") 
            genome << cycle_to_chromosome(chrome)
        end
        return genome
    end

    def two_break_on_genome(p, i, id, j, jd)
        <<-DOC
        2-BreakOnGenome(P, i, i′, j, j′)
             GenomeGraph ← BlackEdges(P) and ColoredEdges(P)
             GenomeGraph ← 2-BreakOnGenomeGraph(GenomeGraph, i, i′, j, j′)
             P ← GraphToGenome(GenomeGraph)
             return P

        Sample Input:
        (+1 -2 -4 +3)
        1, 6, 3, 8
        Sample Output:
        (+1 -2)(-3 +4)    

        DOC

        colored_edges = colored_edges(p)
        colored_edges.each {|edge| print edge.join(",") + "  " }
        puts ""
        mod_colored_edges = two_break_on_genome_graph(colored_edges, i, id, j, jd)
        mod_colored_edges.each {|edge| print edge.join(",") + "  " }
        puts ""
        genome = graph_to_genome(mod_colored_edges)
        # puts genome.join(" ")
        return genome
    end # end of two_break_on_genome function


    def two_break_on_genome_graph(genome_graph, i, id, j, jd)
        <<-DOC
        2-BreakOnGenomeGraph(GenomeGraph, i, i′, j, j′)
             remove colored edges (i, i') and (j, j′) from GenomeGraph
             add colored edges (i, j) and (i′, j') to GenomeGraph
             return GenomeGraph        
        
        Sample Input:
        (2, 4), (3, 8), (7, 5), (6, 1)
        1, 6, 3, 8
        Sample Output:
        (2, 4), (3, 1), (7, 5), (6, 8)
        DOC

        genome_graph.each do |edge|
            if (edge[0] == id && edge[1] == i)
                edge[1] = jd
            elsif (edge[0] == i && edge[1] == id)
                edge[1] = j
            elsif edge[0] == jd && edge[1] == j
                edge[1] = id
            elsif edge[0] == j && edge[1] == jd
                edge[1] = i
            end
        end

        return genome_graph
    end # end of two_break_on_genome_graph function

    def two_break_distance(p, q)
        p_nodes = colored_edges(p)
        q_nodes = colored_edges(q)

        all_nodes = p_nodes + q_nodes

        puts all_nodes.length
        # puts all_nodes.sort.join(" ")

        cycle_count = 0

        loop do 

            break if all_nodes.empty?
            cycle = {}
            # cycle = all_nodes[0]
            cycle[all_nodes[0][0]] = 1
            cycle[all_nodes[0][1]] = 1

            found_a_node = 1

            loop do

                break if found_a_node == 0

                found_a_node = 0
                # all_nodes.dup.each do |node|
                index = 0
                loop do
                    break if all_nodes.length == index
                    node = all_nodes[index]
                    # if (cycle.include?(node[0]) || cycle.include?(node[1]))
                    if (cycle.has_key?(node[0]) || cycle.has_key?(node[1]))
                        # puts all_nodes.length
                        # print node.join(" ") + " "
                        # cycle += node
                        cycle[node[0]] = 1
                        cycle[node[1]] = 1
                        # cycle.uniq!
                        all_nodes.delete(node)
                        found_a_node = 1
                    else
                        index += 1
                    end
                end

            end
            
            # puts cycle.join(" ")
            # puts cycle.keys.length
            cycle_count += 1 unless cycle.empty?
            # puts all_nodes.join(" ")
            # break
        end

        # puts cycle_count
        no_of_blocks = p_nodes.length

        return (no_of_blocks - cycle_count)

    end # end of two_break_distance function


    def colored_edges(p)
        <<-DOC
        ColoredEdges(P)
             Edges ← an empty set
             for each chromosome Chromosome in P
                  Nodes ← ChromosomeToCycle(Chromosome)
                  for j ← 1 to |Chromosome|
                       add the edge (Nodes2j, Nodes2j +1) to Edges
             return Edges        
        The following algorithm constructs ColoredEdges(P) for a genome P. 
        In this pseudocode, we will assume that an n-element array (a1, . . . , an) has an invisible (n + 1)-th element 
        that is equal to its first element, i.e., an+1 = a1.
        DOC

        # puts p
        edges = []
        p.each do |chrome|
            # puts chrome.join(" ")
            nodes = chromosome_to_cyle(chrome)
            # In this pseudocode, we will assume that an n-element array (a1, . . . , an) has an invisible (n + 1)-th element 
            # that is equal to its first element, i.e., an+1 = a1.            
            # So lets add the first element at the end
            nodes.push(nodes.first)
            # puts nodes.join(" ")
            (1..(chrome.length)).each do |j|
                val_2j = 2*j - 1                
                val_2j_1 = (2*j + 1) -1 
                edges << [nodes[val_2j], nodes[val_2j_1]]
            end
            # puts edges.each {|i| puts i.join(",")}
        end
        return edges
    end # end of colored_edges function

    def chromosome_to_cyle(chromosome)
        <<-DOC
        ChromosomeToCycle(Chromosome)
             for j ← 1 to |Chromosome|
                  i ← Chromosomej
                  if i > 0
                       Node(2j−1) ←2i−1
                       Node(2j) ← 2i
                  else
                       Node(2j−1) ← -2i
                       Node(2j) ← -2i−1
             return Nodes        
        DOC

        node = []
        chromosome.each_index do |j|
            i = chromosome[j]
            j += 1 # For j = 0 we will have -ve index for node[]
            if i > 0
                node[2*j-1] = 2*i -1
                node[2*j] = 2*i
            else
                node[2*j-1] = -2*i
                node[2*j] = -2*i - 1
            end
        end
        node.shift # Lets remove the 0'th element in node array as it will always contain a null
        # Because we are starting at j = 1 and we we will always start filling up nodes 2*1 -1 and 2*1 and we never would've 
        # used index 0 in node
        return node
    end # end of chromosome_to_cycle function


    def cycle_to_chromosome(nodes)
        <<-DOC
        CycleToChromosome(Nodes)
             for j ← 1 to |Nodes|/2
                  if Node2j−1 < Node2j
                       Chromosomej ← Node[2j] /2
                  else
                       Chromosomej ← −Node[2j−1]/2
             return Chromosome
        DOC
        chromosome = []

        (1..(nodes.length/2)).each do |j|
            val_2j_1 = 2*j-1 -1 
            val_2j = 2*j -1
            # puts val_2j_1
            # puts val_2j
            if nodes[val_2j_1] < nodes[val_2j]
                chromosome[j-1] = nodes[val_2j]/2
            else
                chromosome[j-1] = -nodes[val_2j_1]/2
            end
        end

        return chromosome
    end # end of cycle_to_chromosome function


    def breakpoints(p)
        <<-DOC
        That consecutive elements (pi pi+1) in permutation P = (p1 ... pn) form an adjacency if pi+1 − pi is equal to 1. 
        By definition, for any positive element k < n, both (k k + 1) and (−(k + 1) −k) are adjacencies. 
        If pi+1 − pi is not equal to 1, then we say that (pi pi+1) is a breakpoint.
        DOC
        no_of_bp = 0

        (0..(p.length-2)).each do |i|
            no_of_bp += 1 if ((p[i+1] - p[i]) != 1)
        end

        return no_of_bp
    end # end of breakpoints function

    def greedy_sorting(p)
        <<-DOC
        GREEDYSORTING(P)
            approxReversalDistance ← 0
            for k = 1 to |P|
                if element k is not sorted
                    apply the k-sorting reversal to P
                    approxReversalDistance ← approxReversalDistance + 1
                if k-th element of P is −k
                    apply the k-sorting reversal to P
                    approxReversalDistance ← approxReversalDistance + 1
            return approxReversalDistance

        Sample Input:
             (-3 +4 +1 +5 -2)

        Sample Output:
             (-1 -4 +3 +5 -2)
             (+1 -4 +3 +5 -2)
             (+1 +2 -5 -3 +4)
             (+1 +2 +3 +5 +4)
             (+1 +2 +3 -4 -5)
             (+1 +2 +3 +4 -5)
             (+1 +2 +3 +4 +5)            
        DOC
        p_a = []
        # puts p.join(" ")
        i = 0
        identity = []
        p.each_index {|i| identity << (i+1)}
        # puts identity.join(" ")
        loop do
            break if p == identity
            new_line = []
            # Lets get the index of the element (i+1) and reverse the array portion from there
            unless (p[i] == (i+1) || p[i] == -(i+1))
                # get the index where our (i+1) element is
                elem_idx = p.index(i+1)
                # our element could have a -ve sign so if elem_idx is nil lets try the -(i+1)
                elem_idx = p.index(-(i+1)) unless elem_idx
                # puts "i #{i} elem_idx #{elem_idx}"
                new_line = p[0..(i-1)] if (i > 0)
                new_line += p[i..elem_idx].reverse.map {|i| i *= -1}

                new_line += p[(elem_idx+1)..(p.length-1)] if (elem_idx < (p.length-1))

                p = new_line
                # puts p.join(" ")
                p_a << p.dup
                # puts p_a.join(" -> ")
            end
            if (p[i] < 0) 
                p[i] *= -1
                # puts p.join(" ")            
                p_a << p.dup
            end
            i += 1
        end
        return p_a
    end # end of greedy_sorting function



    def edit_distance(v,w)
        <<-DOC
        In 1966, Vladimir Levenshtein introduced the notion of the edit distance between two strings as the minimum number of edit operations 
        needed to transform one string into another.
        Sample Input:
        PLEASANTLY
        MEANLY
        Sample Output:
        5

        # http://en.wikipedia.org/wiki/Levenshtein_distance
        DOC


        s = []
        (0..(v.length)).each {|i| s[i] = []}
        (0..(v.length)).each {|i| s[i][0] = i}
        (0..(w.length)).each {|j| s[0][j] = j}

        (1..(v.length)).each do |i|
            (1..(w.length)).each do |j|
                if v[i] == w[j] 
                    s[i][j] = s[i-1][j-1]
                else
                    s[i][j] = [s[i-1][j]+1, s[i][j-1]+1, s[i-1][j-1]+1].min
                end
            end
        end
        puts s[v.length][w.length]
    end # end of edit_distance function


    def dag_longpath(source, sink, weight_edge_graph, graph)
        <<-DOC
        The backtracking method can be generalized to construct a longest path in any DAG. Whenever we compute S(b) as
        S(b) = max (all predecessors a of node b) {S(a) + weight_of_edge_from_a_to_b }
        we simply need to store a predecessor of b that was used in the computation of sb so that we can backtrack later on. 
        You are now ready to use backtracking to find the longest path in an arbitrary DAG.
        DOC
        topo_graph = Marshal.load(Marshal.dump(graph))
        topo_nodes = topological_ordering(topo_graph)
        # puts topo_nodes.join(" ")

        puts graph
        keys = graph.keys
        values = graph.values.flatten.uniq
        # puts keys.join(" ")
        # puts values.join(" ")
        nodes_with_indegree_zero = keys - values
        puts nodes_with_indegree_zero.join(" ")

        puts weight_edge_graph.join("  ")
        nodes_with_indegree_zero.each do |node|
            unless node == source
                topo_nodes.delete(node) 
                weight_edge_graph.delete_if {|entry| /^#{node}->/.match(entry)}
            end
        end

        # nodes_to_consider = [source]
        # next_node = source
        # index = 1
        # loop do
        #     break if next_node == sink
        #     if graph[next_node]
        #         graph[next_node].each do |node|
        #             # puts node
        #             nodes_to_consider << node unless nodes_to_consider.include?(node)
        #         end
        #     end
        #     next_node = nodes_to_consider[index]
        #     index += 1
        # end
        # puts "nodes_to_consider: " + nodes_to_consider.join(" ")
        # topo_nodes.keep_if {|node| nodes_to_consider.include?(node)}
        puts "topo nodes: " + topo_nodes.join(" ")

        puts weight_edge_graph.join("  ")
        incoming_nodes = {}
        incoming_nodes[topo_nodes[0]] = [{"pred_node" => nil, "weight" => 0}]


        best_pred = {}
        topo_nodes.each do |node|
            # incoming_nodes[node] = [{"pred_node" => nil, "weight" => 0}]
            best_pred[node] = {"pred_node" => nil, "weight" => 0}
        end

        # puts weight_edge_graph.join("  ")
        (1..(topo_nodes.length-1)).each do |i|
            node = topo_nodes[i]
            puts "Node: " + node
            all_weights = []
            weight_edge_graph.each do |entry| 
                # puts "Entry: " + entry
                # puts best_pred
                if /(\d+)->#{node}:(\d+)/.match(entry) 
                    incoming_nodes[node] = [] unless incoming_nodes[node]
                    pred_node_matched = /(\d+)->#{node}:(\d+)/.match(entry)[1]
                    pred_node_weight = /(\d+)->#{node}:(\d+)/.match(entry)[2].to_i
                    # puts /(\d+)->#{node}:(\d+)/.match(entry)[1]
                    # puts /(\d+)->#{node}:(\d+)/.match(entry)[2]
                    incoming_nodes[node] << {"pred_node" => pred_node_matched, "weight" => pred_node_weight}
                    all_weights << (pred_node_weight + best_pred[pred_node_matched]["weight"])
                end
            end
            # puts all_weights.join(" ")
            if incoming_nodes[node] 
                # puts incoming_nodes[node]
                incoming_nodes[node].each do |pred_node_entry|
                    # puts all_weights.join(" ")
                    if ((best_pred[pred_node_entry["pred_node"]]["weight"] + pred_node_entry["weight"]) == all_weights.max)
                        best_pred[node]["weight"] = all_weights.max
                        best_pred[node]["pred_node"] = pred_node_entry["pred_node"]
                    end
                end
            end
        end
        # puts incoming_nodes
        # puts best_pred
        length = best_pred[sink]["weight"]
        path = [sink]
        node = sink
        loop do
            puts node
            break if node == source
            node = best_pred[node]["pred_node"]
            path << node
        end
        return length, path.reverse
    end # end of dag_longpath function


    def local_alignment(v, w, sigma)
        val, i, j, backtrack = lcs_backtrack_pam(v, w, sigma)
        backtrack.each do |bk|
            # puts bk.join(" ")
        end
        puts val
        output_lcs_pam(backtrack, v, w, i, j)
    end # end of local_alignment function


    def output_lcs_pam(backtrack, v, w, i, j)

        v_str = ""
        w_str = ""

        loop do
            # break if (i == 0 || j == 0)
            break if (i == 0)
            # puts i.to_s + " " + j.to_s + " " + backtrack[i][j]
            if backtrack[i][j] == "free"
                break
            elsif backtrack[i][j] == "down"
                v_str += v[i-1]
                w_str += "-"
                i -= 1
            elsif backtrack[i][j] == "right"
                v_str += "-"
                w_str += w[j-1]
                j -= 1
            else
                v_str += v[i-1]
                w_str += w[j-1]
                i -= 1
                j -= 1
            end
        end

        puts v_str.reverse
        puts w_str.reverse

    end # end of output_lcs_pam function

    def lcs_backtrack_pam(v, w, sigma)
        <<-DOC
        LCSBACKTRACK(v, w)
        for i ← 0 to |v|
            si, 0 ← 0
        for j ← 0 to |w| 
            s0, j ← 0
        for i ← 1 to |v|
            for j ← 1 to |w|
                si, j ← max{si-1, j, si,j-1, si-1, j-1 + 1 (if vi = wj)}
                if si,j = si-1,j
                    Backtracki, j ← "↓"
                if si, j = si, j-1
                    Backtracki, j ← "→"
                if si, j = si-1, j-1 + 1
                    Backtracki, j ← "↘"
        return Backtrack
        DOC
        # puts v
        # puts w

        pam = pam_250()
        s = []
        backtrack = []

        max_score = 0
        max_i = 0
        max_j = 0

        (0..(v.length)).each {|i| s[i] = []}
        (0..(v.length)).each {|i| backtrack[i] = []}

        (0..(v.length)).each do |i| 
            # s[i][0] = i*(-sigma)
            s[i][0] = 0
            backtrack[i][0] = "down"
        end

        (0..(w.length)).each do |j| 
            # s[0][j] = j*(-sigma)
            s[0][j] = 0
            backtrack[0][j] = "right"
        end

        (1..(v.length)).each do |i|
            (1..(w.length)).each do |j|
                # print  "i #{i} j #{j} #{v[i-1]} #{w[j-1]} #{pam[v[i-1]][w[j-1]]} "
                down_val = s[i-1][j] - sigma
                right_val = s[i][j-1] - sigma
                diag_val = s[i-1][j-1] + pam[v[i-1]][w[j-1]]
                s[i][j] = [down_val, right_val, diag_val].max
                if s[i][j] <= 0
                    s[i][j] = 0
                    backtrack[i][j] = "free"
                elsif s[i][j] == down_val
                    backtrack[i][j] = "down"
                elsif s[i][j] == right_val
                    backtrack[i][j] = "right"
                elsif s[i][j] == diag_val
                    backtrack[i][j] = "diag"
                end
                # print [down_val, right_val, diag_val].join(" ") + "  "
                # print  "S: " + s[i][j].to_s + "  "
                # print "Backtrack: " + backtrack[i][j] + " \n"
                if max_score <= s[i][j]
                    max_score = s[i][j]
                    max_i = i
                    max_j = j
                end
            end            
        end

        return [s[max_i][max_j], max_i, max_j, backtrack]
    end # end of lcs_backtrack_pam function   


    def global_alignment(v, w, sigma)
        # https://class.coursera.org/bioinformatics-002/forum/thread?thread_id=634
        # https://drive.google.com/file/d/0B-C9of9KPRGKa09tX2pkb0t0TEU/view
        val, backtrack = lcs_backtrack_blosum(v, w, sigma)
        backtrack.each do |bk|
            # puts bk.join(" ")
        end
        puts val
        out_v, out_w = output_lcs_blosum(backtrack, v, w, v.length, w.length)
        return [out_v, out_w]
    end # end of global_alignment function


    def output_lcs_blosum(backtrack, v, w, i, j)

        v_str = ""
        w_str = ""

        loop do
            # break if (i == 0 || j == 0)
            break if (i == 0)
            # puts i.to_s + " " + j.to_s + " " + backtrack[i][j]
            if backtrack[i][j] == "down"
                v_str += v[i-1]
                w_str += "-"
                i -= 1
            elsif backtrack[i][j] == "right"
                v_str += "-"
                w_str += w[j-1]
                j -= 1
            else
                v_str += v[i-1]
                w_str += w[j-1]
                i -= 1
                j -= 1
            end
        end

        # puts v_str.reverse
        # puts w_str.reverse
        return [v_str.reverse, w_str.reverse]
        # return [v_str, w_str]
        <<-DOC
         OUTPUTLCS(backtrack, v, i, j)
                if i = 0 or j = 0
                    return
                if backtracki, j = "↓"
                    OUTPUTLCS(backtrack, v, i - 1, j)
                else if backtracki, j = "→"
                    OUTPUTLCS(backtrack, v, i, j - 1)
                else
                    OUTPUTLCS(backtrack, v, i - 1, j - 1)
                    output vi        
        DOC

        # return if (i == 0 || j == 0)
        # if backtrack[i][j] == "del"
        #     output_lcs_blosum(backtrack, v, w, i-1, j) 
        #     print "-"
        # elsif backtrack[i][j] == "add"
        #     output_lcs_blosum(backtrack, v, w, i, j-1) 
        #     print "-"
        # # elsif backtrack[i][j] == "mismatch"
        # #     output_lcs(backtrack, v, w, i-1, j-1) 
        # #     print "-"
        # else
        #     output_lcs_blosum(backtrack, v, w, i-1, j-1) 
        #     print v[i-1]
        #     print w[j-1]
        # end
    end # end of output_lcs_blosum function

    def lcs_backtrack_blosum(v, w, sigma)
        <<-DOC
        LCSBACKTRACK(v, w)
        for i ← 0 to |v|
            si, 0 ← 0
        for j ← 0 to |w| 
            s0, j ← 0
        for i ← 1 to |v|
            for j ← 1 to |w|
                si, j ← max{si-1, j, si,j-1, si-1, j-1 + 1 (if vi = wj)}
                if si,j = si-1,j
                    Backtracki, j ← "↓"
                if si, j = si, j-1
                    Backtracki, j ← "→"
                if si, j = si-1, j-1 + 1
                    Backtracki, j ← "↘"
        return Backtrack
        DOC
        # puts v
        # puts w

        blosum = blosum_62()
        s = []
        backtrack = []
        (0..(v.length)).each {|i| s[i] = []}
        (0..(v.length)).each {|i| backtrack[i] = []}

        # (0..(v.length-1)).each {|i| s[i][0] = 0}
        (0..(v.length)).each do |i| 
            s[i][0] = i*(-sigma)
            backtrack[i][0] = "down"
        end
        # (0..(w.length-1)).each {|j| s[0][j] = 0}
        (0..(w.length)).each do |j| 
            s[0][j] = j*(-sigma)
            backtrack[0][j] = "right"
        end

        (1..(v.length)).each do |i|
            (1..(w.length)).each do |j|
                # puts "i #{i} j #{j} #{v[i-1]} #{w[j-1]} #{blosum[v[i-1]][w[j-1]]}"
                down_val = s[i-1][j] - sigma
                right_val = s[i][j-1] - sigma
                diag_val = s[i-1][j-1] + blosum[v[i-1]][w[j-1]]
                s[i][j] = [down_val, right_val, diag_val].max
                if s[i][j] == down_val
                    backtrack[i][j] = "down"
                elsif s[i][j] == right_val
                    backtrack[i][j] = "right"
                elsif s[i][j] == diag_val
                    backtrack[i][j] = "diag"
                end

            end            
        end

        return [s[v.length][w.length], backtrack]
    end # end of lcs_backtrack_blosom function    

    def output_lcs(backtrack, v, i, j)
        <<-DOC
         OUTPUTLCS(backtrack, v, i, j)
                if i = 0 or j = 0
                    return
                if backtracki, j = "↓"
                    OUTPUTLCS(backtrack, v, i - 1, j)
                else if backtracki, j = "→"
                    OUTPUTLCS(backtrack, v, i, j - 1)
                else
                    OUTPUTLCS(backtrack, v, i - 1, j - 1)
                    output vi        
        DOC

        return if (i == 0 || j == 0)
        if backtrack[i][j] == "del"
            output_lcs(backtrack, v, i-1, j) 
            # print "-"
        elsif backtrack[i][j] == "add"
            output_lcs(backtrack, v, i, j-1) 
            # print "-"
        # elsif backtrack[i][j] == "mismatch"
        #     output_lcs(backtrack, v, i-1, j-1) 
        #     print "-"
        else
            output_lcs(backtrack, v, i-1, j-1) 
            print v[i-1]
        end
    end

    def lcs_backtrack(v,w)
        <<-DOC
        LCSBACKTRACK(v, w)
        for i ← 0 to |v|
            si, 0 ← 0
        for j ← 0 to |w| 
            s0, j ← 0
        for i ← 1 to |v|
            for j ← 1 to |w|
                si, j ← max{si-1, j, si,j-1, si-1, j-1 + 1 (if vi = wj)}
                if si,j = si-1,j
                    Backtracki, j ← "↓"
                if si, j = si, j-1
                    Backtracki, j ← "→"
                if si, j = si-1, j-1 + 1
                    Backtracki, j ← "↘"
        return Backtrack
        DOC
        s = []
        backtrack = []
        # (0..(v.length-1)).each {|i| s[i] = []}
        (0..(v.length)).each {|i| s[i] = []}
        # (0..(v.length-1)).each {|i| backtrack[i] = []}
        (0..(v.length)).each {|i| backtrack[i] = []}

        # (0..(v.length-1)).each {|i| s[i][0] = 0}
        (0..(v.length)).each {|i| s[i][0] = 0}
        # (0..(w.length-1)).each {|j| s[0][j] = 0}
        (0..(w.length)).each {|j| s[0][j] = 0}

        # (0..(v.length-1)).each do |i|
        (1..(v.length)).each do |i|
            # (0..(w.length-1)).each do |j|
            (1..(w.length)).each do |j|
                # puts "i #{i} j #{j} #{v[i-1]} #{w[j-1]}"
                val = 0
                val = s[i-1][j-1] + 1 if (v[i-1] == w[j-1])
                s[i][j] = [s[i-1][j], s[i][j-1], val].max
                if s[i][j] == s[i-1][j]
                    backtrack[i][j] = "del" 
                elsif s[i][j] == s[i][j-1]
                    backtrack[i][j] = "add" 
                elsif (s[i][j] == val)
                    backtrack[i][j] = "match" 
                end
            end            
        end
        puts s[v.length][w.length]
        return backtrack
    end # end of lcs_backtrack function

    def topological_ordering(graph)
        <<-DOC
        TOPOLOGICALORDERING(Graph)
            List ← empty list
            Candidates ← set of all nodes in Graph with no incoming edges
            while Candidates is non-empty
                select an arbitrary node a from Candidates
                add a to the end of List and remove it from Candidates
                for each outgoing edge from a to another node b
                    remove edge (a, b) from Graph
                    if b has no other incoming edges 
                        add b to Candidates
            if Graph has edges that have not been removed
                return "the input graph is not a DAG"
            else return List        
        DOC

        list = []
        keys = graph.keys
        values = graph.values.flatten.uniq
        candidates = keys - values
        loop do
            break if candidates.empty?
            a = candidates[0]
            list << a
            candidates.delete(a)
            loop do
                break unless graph[a]
                b = graph[a][0]
                graph[a].delete(b)
                graph.delete(a) if graph[a].empty?
                candidates << b unless graph.values.flatten.include?(b)
            end
        end

        unless graph.empty?
            return "Graph is not a DAG"            
        end

        return list
    end # end of topological_ordering function

    def manhattan_tourist(n, m, down, right)
        <<-DOC
        We now have the outline of a dynamic programming algorithm for finding the length of a longest path in the Manhattan Tourist Problem, 
        called MANHATTANTOURIST. In the following pseudocode, downi, j and righti, j are the respective weights of the vertical and 
        horizontal edges entering node (i, j). We denote the matrices holding (downi, j) and (righti, j) as Down and Right, respectively.

        Input: Integers n and m, followed by an n × (m + 1) matrix Down and an (n + 1) × m matrix Right.
        The two matrices are separated by the - symbol.
        Output: The length of a longest path from source (0, 0) to sink (n, m) in the n × m rectangular grid
        whose edges are defined by the matrices Down and Right.

        MANHATTANTOURIST(n, m, Down, Right)
        s0, 0 ← 0
        for i ← 1 to n
            s(i, 0) ← s(i-1, 0) + down(i, 0)
        for j ← 1 to m
            s(0, j) ← s(0, j−1) + right(0, j)
        for i ← 1 to n
            for j ← 1 to m
                s(i, j) ← max{s(i - 1, j) + down(i, j), s(i, j - 1) + right(i, j)}
        return s(n, m)
        DOC

        s = []
        (0..(n)).each {|i| s[i] = []}
        s[0][0] = 0        
        (1..(n)).each {|i| s[i][0] = s[i-1][0] + down[i-1][0] }
        (1..(m)).each {|j| s[0][j] = s[0][j-1] + right[0][j-1] }
        # puts s.join(" ")
        (1..n).each do |i|
            (1..m).each do |j|
                s[i][j] = [s[i-1][j] + down[i-1][j] , s[i][j-1] + right[i][j-1] ].max
            end
        end
        return s[n][m]
    end # end of manhattan_tourist function

    def dp_change(money, coins) 
        <<-DOC
        Input: An integer money and an array Coins = (coin1, ..., coind).
        Output: The minimum number of coins with denominations Coins that changes money.
        
        DPCHANGE(money, Coins)
            MinNumCoins(0) ← 0
            for m ← 1 to money
                MinNumCoins(m) ← ∞
                    for i ← 1 to |Coins|
                        if m ≥ coini
                            if MinNumCoins(m - coini) + 1 < MinNumCoins(m)
                                MinNumCoins(m) ← MinNumCoins(m - coini) + 1
            output MinNumCoins(money)        
        
        DOC

        # Since we can't assign infinity, lets assign the max possible coins + 1 as inifinity
        # Get the least value coin and divide our money with that 
        infinity = (money/coins[-1]) + 1
        min_num_coins = {}
        min_num_coins[0] = 0
        (1..money).each do |m|
            min_num_coins[m] = infinity
            coins.each do |coin|
                if (m >= coin)
                    if (min_num_coins[m-coin] + 1) < min_num_coins[m]
                        min_num_coins[m] = min_num_coins[m-coin] + 1
                    end
                end
            end
        end
        return min_num_coins[money]
    end # end of dp_change function

    def universal_string(k)
        bin_kmers = []
        (0..(2**k-1)).each {|n| bin_kmers << ( "%0#{k}b" % n)}
        puts bin_kmers.join(" ")
        graph = debruijn_graph_from_kmers(bin_kmers)
        puts graph.to_s
        path = eulerian_path(graph)
        # puts path.join(" ")

        # If there is consecutive same values remove those
        # 000 000 001 010 100 001 011 110 101 010 101 011 111 111 110 100 000

        return genome_from_path(path)
    end # end of universal_string function


    def eulerian_path(graph)
        # We can just use the eulerian cycle problem by just choosing the 
        # first node which has greater outdegree than indegree

        node = graph.keys[0]
        in_degrees = {}
        graph.values.flatten.each do |val|
            in_degrees[val] = 0 unless in_degrees.has_key?(val)
            in_degrees[val] += 1
        end
        # puts in_degrees.to_s
        graph.each do |k,v|
            if in_degrees.has_key?(k) and v.length > in_degrees[k]
                node = k
                break
            end
            unless in_degrees.has_key?(k)
                node = k
                break
            end
        end
        stack = [node]
        cycle = []
        while (not stack.empty?)
            node = stack[-1]
            if graph.has_key?(node)
                new_node = graph[node][0]
                # puts "New node:" + new_node
                stack << new_node
                if graph[node].length > 1
                    graph[node].delete(new_node)
                else
                    graph.delete(node)
                end
            else
                cycle.unshift(stack.pop)
                # puts stack.join(" ")
            end
        end

        return cycle     
    end # end of eulerian_path function


    def paired_composition(text, k, d)
        <<-DOC
        Given a string Text, a (k,d)-mer is a pair of k-mers in Text separated by distance d. 
        We use the notation (Pattern1|Pattern2) to refer to a a (k,d)-mer whose k-mers are Pattern1 and Pattern2. 
        For example, (ATG|GGG) is a (3,4)-mer in TAATGCCATGGGATGTT. The (k,d)-mer composition of Text, denoted PairedCompositionk,d(Text), 
        is the collection of all (k,d)- mers in Text (including repeated (k,d)-mers). For example, 
        here is PairedComposition3,1(TAATGCCATGGGATGTT):

            TAA GCC          
             AAT CCA         
              ATG CAT        
               TGC ATG       
                GCC TGG      
                 CCA GGG     
                  CAT GGA    
                   ATG GAT   
                    TGG ATG  
                     GGG TGT 
                      GGA GTT
            TAATGCCATGGGATGTT
            Since the order of (3,1)-mers is unknown, we list them according to the lexicographic order of the 6-mers formed by their concatenated 3-mers:

            (AAT|CCA) (ATG|CAT) (ATG|GAT) (CAT|GGA) (CCA|GGG) (GCC|TGG) (GGA|GTT) (GGG|TGT) (TAA|GCC) (TGC|ATG) (TGG|ATG)
        DOC

        strings = []
        (0..(text.length - 2*k - d)).each do |i|
            strings << text.slice(i,2*k+d)
        end
        strings.sort!
        # AATGCCA ATGCCAT ATGGGAT CATGGGA CCATGGG GCCATGG GGATGTT GGGATGT TAATGCC TGCCATG TGGGATG

        # Now lets replace the 'd' characters in the middle with a ',' so we can easily parse the output to represent the paried composition
        strings.map! do |str| 
            str.slice!(k, d)
            str.insert(k, '|')
        end

        return strings
    end # end of paired_composition function   


    def eulerian_cycle(graph)
        <<-DOC
        EULERIANCYCLE(Graph)
        form a cycle Cycle by randomly walking in Graph (don't visit the same edge twice!)
        while there are unexplored edges in Graph
            select a node newStart in Cycle with still unexplored edges
            form Cycle’ by traversing Cycle (starting at newStart) and then randomly walking 
            Cycle ← Cycle’
        return Cycle

        http://www.ms.uky.edu/~lee/ma515fa10/euler.pdf

        This is an algorithm to find an Eulerian circuit in a connected graph in which every vertex has even degree.
            1. Choose any vertex v and push it onto a stack. Initially all edges are unmarked.
            2. While the stack is nonempty, look at the top vertex, u, on the stack. 
                If u has an unmarked incident edge, say, to a vertex w, then push w onto the stack 
                and mark the edge uw. On the other hand, if u has no unmarked incident edge, then pop u off the stack and print it.
            
            When the stack is empty, you will have printed a sequence of vertices that correspond to an Eulerian circuit.

        DOC

        node = graph.keys[0]
        stack = [node]
        cycle = []
        while (not stack.empty?)
            node = stack[-1]
            if graph.has_key?(node)
                new_node = graph[node][0]
                # puts "New node:" + new_node
                stack << new_node
                if graph[node].length > 1
                    graph[node].delete(new_node)
                else
                    graph.delete(node)
                end
            else
                cycle.unshift(stack.pop)
                # puts stack.join(" ")
            end
        end

        # return cycle.join("->")
        return cycle

    end # end of eulerian_cycle function

    def string_reconstruction(k, kmers)
        graph = debruijn_graph_from_kmers(kmers)
        # puts graph.to_s
        path = eulerian_path(graph)
        # puts path.join(" ")
        return genome_from_path(path)
    end # end of string_reconstruction function

    def string_spelled_by_gapped_patterns(patterns, k, d)
        <<-DOC
        StringSpelledByGappedPatterns(GappedPatterns, k, d)
            FirstPatterns ← the collection of initial k-mers from GappedPatterns
            SecondPatterns ← the collection of terminal k-mers from GappedPatterns
            PrefixString ← StringSpelledByPatterns(FirstPatterns, k)
            SuffixString ← StringSpelledByPatterns(SecondPatterns, k)
            for i = k + d + 1 to |PrefixString|
                if the i-th symbol in PrefixString does not equal the (i - k - d)-th symbol in SuffixString
                    return "there is no string spelled by the gapped patterns"
            return PrefixString concatenated with the last k + d symbols of SuffixString
        DOC
        first_patterns = []
        second_patterns = []
        patterns.each do |pat|
            first_patterns << pat.split("|")[0]
            second_patterns << pat.split("|")[1]
        end
        # puts first_patterns
        # puts second_patterns
        prefix_string = genome_from_path(first_patterns)
        # puts prefix_string
        suffix_string = genome_from_path(second_patterns)
        # puts suffix_string
        ((k+d+1)..(prefix_string.length-1)).each do |i|
            return nil unless prefix_string[i] == suffix_string[i-k-d]
        end
        return (prefix_string + suffix_string.slice(suffix_string.length-k-d, k+d))
    end # end of string_spelled_by_gapped_patterns function

    def string_reconstruction_from_read_pairs(k, d, kmers)
        graph = debruijn_graph_from_paired_kmers(kmers, k, d)
        # puts graph.to_s
        path = eulerian_path(graph)
        # puts path.join(" ")
        return string_spelled_by_gapped_patterns(path, k, d)
    end # end of string_reconstruction function


    def path_rearrange(path,node)
        index = path.index(node)
        # If the index is 0, the cycle is going to be same so just return path
        return path if index == 0
        # Path 0 3 2 1 0
        # Node 2
        # Arranged path = 2 1 0 3 2
        new_path = path.slice(index, path.length-index)
        new_path.push(path.slice(1,index))
        return new_path
    end

    def string_composition(k, text)
        <<-DOC
        String Composition Problem: Generate the k-mer composition of a string.
             Input: An integer k and a string Text.
             Output: Compositionk(Text), where the k-mers are arranged in lexicographic order.
        DOC

        return kmers(text,k).sort
    end # end of string_composition function

    def genome_from_path(kmers_a)
        <<-DOC
        String Spelled by a Genome Path Problem. Reconstruct a string from its genome path.
         Input: A sequence of k-mers Pattern1, … ,Patternn such that the last k - 1 symbols of Patterni are
                    equal to the first k-1 symbols of Patterni+1 for 1 ≤ i ≤ n-1.
         Output: A string Text of length k+n-1 such that the i-th k-mer in Text is equal to Patterni  (for 1 ≤ i ≤ n).        
            
          Input: ["TAA", "AAT", "ATG"]
          Output: TAATG
        DOC
        # If we assign the first kmer to the genome and keep on adding the last character from all the next kmers we
        # should get the string
        genome = kmers_a.shift # shift will remove first element from array and returns it
        kmers_a.each do |kmer|
            genome += kmer[-1]
        end
        return genome
    end # end of genome_from_path function

    def debruijn_graph(k, text)
        <<-DOC
        In general, given a genome Text, PathGraphk(Text) is the path consisting of |Text| - k + 1 edges, 
        where the i-th edge of this path is labeled by the i-th k-mer in Text and the i-th node of the path 
        is labeled by the i-th (k - 1)-mer in Text. 
        The de Bruijn graph DeBruijnk(Text) is formed by gluing identically labeled nodes in PathGraphk(Text).

        De Bruijn Graph from a String Problem: Construct the de Bruijn graph of a string.
             Input: An integer k and a string Text.
             Output: DeBruijnk(Text).        
        DOC
        kmers_a = kmers(text, k-1)

        graph = {}
        (0..(kmers_a.length-1-1)).each do |i|
            kmer = kmers_a[i]
            graph[kmer] = [] unless graph.has_key?(kmer)
            graph[kmer] << kmers_a[i+1] 
        end
        # puts graph
        return graph
    end # end of debgruijn_graph function

    def maximal_non_branching_paths(graph)
        <<-DOC
        A node v in a directed graph Graph is called a 1-in-1-out node if its indegree and outdegree are both equal to 1, 
        i.e., in(v) = out(v) = 1.  We can rephrase the definition of a "maximal non-branching path" from the main text 
        as a path whose internal nodes are 1-in-1-out nodes and whose initial and final nodes are not 1-in-1-out nodes.  
        Also, note that the definition from the main text does not handle the special case when Graph has a connected 
        component that is an isolated cycle, in which all nodes are 1-in-1-out nodes.

        The MaximalNonBranchingPaths pseudocode below generates all non-branching paths in a graph. 
        It iterates through all nodes of the graph that are not 1-in-1-out nodes and generates all non-branching 
        paths starting at each such node.  In a final step, MaximalNonBranchingPaths finds all isolated cycles in the graph.

        MaximalNonBranchingPaths(Graph)
            Paths ← empty list
            for each node v in Graph
                if v is not a 1-in-1-out node
                    if out(v) > 0
                        for each outgoing edge (v, w) from v
                            NonBranchingPath ← the path consisting of the single edge (v, w)
                            while w is a 1-in-1-out node
                                extend NonBranchingPath by the outgoing edge (w, u) from w 
                                w ← u
                            add NonBranchingPath to the set Paths
            for each isolated cycle Cycle in Graph
                add Cycle to Paths
            return Paths        
        DOC

        paths = []

        in_degrees = {}
        out_degrees = {}
        # Calculating the indegree and outdegree for each node
        graph.values.flatten.each do |val|
            in_degrees[val] = 0 unless in_degrees.has_key?(val)
            in_degrees[val] += 1
        end
        graph.keys.each do |val|
            out_degrees[val] = graph[val].length
        end

        non_branch_path = ""

        one_in_one_out_nodes = []
        in_degrees.each do |k,v|
            if v == 1
                one_in_one_out_nodes << k if (out_degrees[k] == 1)
            end
        end

        # puts one_in_one_out_nodes
        nodes = []
        nodes = graph.keys
        graph.values.flatten.each {|v| nodes << v}
        nodes.uniq!

        v_index = 0
        loop do
            break if v_index == nodes.length
            v = nodes[v_index]

            if (one_in_one_out_nodes.include?(v))
                v_index += 1
                next
            else
                if (out_degrees[v] && (out_degrees[v] > 0))
                    loop do
                        break unless graph[v]
                        break if graph[v].empty?
                        w = graph[v][0]
                        non_branch_path = v
                        non_branch_path += w[-1]
                        graph[v].delete_at(graph[v].index(w))
                        graph.delete(v) if graph[v].length == 0
                        while (w && graph[w] && in_degrees[w] && out_degrees[w] && in_degrees[w] == 1 && out_degrees[w] == 1) do
                            u = graph[w][0]
                            non_branch_path += u[-1]
                            graph.delete(w)
                            w = u
                        end
                        paths << non_branch_path
                        # puts "path:" + non_branch_path
                        non_branch_path = ""
                    end
                else
                    v_index += 1
                    next
                end
            end
            v_index += 1
        end
        # puts paths
        paths.sort
    end # end of maximal_non_branching_paths function

    def contig_generation(kmers_a)
        graph = debruijn_graph_from_kmers(kmers_a)
        # puts graph.keys.length
        # puts graph.values.flatten.length
        # puts graph
        return maximal_non_branching_paths(graph)
    end # end of contig_generation function


    def debruijn_graph_from_kmers(kmers_a)
        <<-DOC
        DeBruijn Graph from k-mers Problem: Construct the de Bruijn graph from a set of k-mers.
             Input: A collection of k-mers Patterns.
             Output: The adjacency list of the de Bruijn graph DeBruijn(Patterns).        
        DOC

        graph = {}
        kmers_a.each do |kmer|
            kmer_p = kmer_prefix(kmer)
            kmer_s = kmer_suffix(kmer)
            graph[kmer_p] = [] unless graph.has_key?(kmer_p)
            graph[kmer_p] << kmer_s
        end

        return graph
    end # end of debruijn_graph_from_kmers function


    def debruijn_graph_from_paired_kmers(kmers_a, k, d)
        <<-DOC
        DeBruijn Graph from k-mers Problem: Construct the de Bruijn graph from a set of k-mers.
             Input: A collection of k-mers Patterns.
             Output: The adjacency list of the de Bruijn graph DeBruijn(Patterns).        
        DOC

        graph = {}
        kmers_a.each do |kmer|
            kmer_p = paired_kmer_prefix(kmer, k, d)
            kmer_s = paired_kmer_suffix(kmer, k, d)
            graph[kmer_p] = [] unless graph.has_key?(kmer_p)
            graph[kmer_p] << kmer_s
        end

        return graph
    end # end of debruijn_graph_from_kmers function


    def overlap_graph(kmer_a)
        <<-DOC
        CODE CHALLENGE: Solve the Overlap Graph Problem (restated below).
            Input: A collection Patterns of k-mers.
            Output: The overlap graph Overlap(Patterns), in the form of an adjacency list.
        DOC
        # puts kmer_a
        graph = {}
        (0..(kmer_a.length-1)).each do |i|
            (0..(kmer_a.length-1)).each do |j|
                unless i == j
                    # puts kmer_a[i] + " " + kmer_a[j]
                    # puts "suffix " + kmer_suffix(kmer_a[i])
                    # puts "prefix " + kmer_prefix(kmer_a[j])
                    if kmer_suffix(kmer_a[i]) == kmer_prefix(kmer_a[j])
                        # puts "matched"
                        if graph.has_key?(kmer_a[i])
                            graph[kmer_a[i]] << kmer_a[j]
                        else
                            graph[kmer_a[i]] = [kmer_a[j]]
                        end
                    end
                end
            end
        end

        # puts graph
        return graph
    end # end of overlap_graph function

    def median_strings(dna, k)
        # MedianString(Dna, k)
        # distance ← ∞
        # for i ←0 to 4k −1
        #     Pattern ← NumberToPattern(i, k)
        #     if distance > DistanceBetweenPatternAndStrings(Pattern, Dna)
        #         distance ← DistanceBetweenPatternAndStrings(Pattern, Dna)
        #         Median ← Pattern
        # return Median        

        pattern = number_to_pattern(0, k)
        median = pattern
        distance = min_distance_pattern_dna(pattern, dna)
        # puts distance
        median_arr = []
        (1..(4**k-1)).each do |i|
            pattern = number_to_pattern(i, k)
            # puts "GAC: " + min_distance_pattern_dna(pattern, dna).to_s if pattern == "GAC"
            if distance >= min_distance_pattern_dna(pattern, dna)
                distance = min_distance_pattern_dna(pattern, dna)
                # puts distance
                # median = pattern
                median_arr << pattern
            end
        end

        return median_arr
    end # End of median_strings function

    def median_string(dna, k)
        # MedianString(Dna, k)
        # distance ← ∞
        # for i ←0 to 4k −1
        #     Pattern ← NumberToPattern(i, k)
        #     if distance > DistanceBetweenPatternAndStrings(Pattern, Dna)
        #         distance ← DistanceBetweenPatternAndStrings(Pattern, Dna)
        #         Median ← Pattern
        # return Median        

        pattern = number_to_pattern(0, k)
        median = pattern
        distance = min_distance_pattern_dna(pattern, dna)
        # puts distance
        (1..(4**k-1)).each do |i|
            pattern = number_to_pattern(i, k)
            # puts "GAC: " + min_distance_pattern_dna(pattern, dna).to_s if pattern == "GAC"
            if distance > min_distance_pattern_dna(pattern, dna)
                distance = min_distance_pattern_dna(pattern, dna)
                # puts distance
                median = pattern
            end
        end

        return median
    end # End of median_string function

    def min_distance_pattern_dna(pattern, dna)
        pattern.upcase!
        # d(Pattern,Dna) = ∑  d(Pattern,Dna[i]).
        min_dist = 0
        dna.each do |text|
            text.upcase!
            min_dist += min_distance_pattern_text(pattern, text)
        end
        return min_dist
    end # End of min_distance_pattern_dna function

    def min_distance_pattern_text(pattern, text)
        pattern.upcase!
        text.upcase!
        # d(Pattern,Text)= min of all k-mers Pattern' in Text calculate the HammingDistance(Pattern,Pattern′)
        min_dist = hamming_distance(pattern, text.slice(0,pattern.length))
        (1..(text.length-pattern.length)).each do |i|
            kmer = text.slice(i,pattern.length)
            min_dist = hamming_distance(pattern, kmer) if min_dist > hamming_distance(pattern, kmer)
        end
        return min_dist
    end # End of min_distance_pattern_text function

    def probability_entropy(prob_arr)
        # Given a probability array get the entropy value
        # For example, the entropy of the probability distribution (0.2, 0.6, 0.0, 0.2) corresponding to the 2nd column of the NF-κB profile matrix is
        # −(0.2log(2)0.2 + 0.6log(2)0.6 + 0.0log(2)0.0 + 0.2log(2)0.2) ~ 1.371        
        entropy = 0.0
        prob_arr.each do |val|
                entropy += (val*Math.log2(val)) unless val == 0.0
        end
        return entropy * -1.0
    end # end of probability_entropy function

    def motifs_entropy(motifs)
        # Entropy is a measure of the uncertainty of a probability distribution (p1, …, pN), and is defined as follows:
        # H(p1,…,pN)= − for i = 1 to N { pi·log2pi }
        # For example, the entropy of the probability distribution (0.2, 0.6, 0.0, 0.2) corresponding to the 2nd column of the NF-κB profile matrix is

        # −(0.2log(2)0.2 + 0.6log(2)0.6 + 0.0log(2)0.0 + 0.2log(2)0.2) ~ 1.371

        # We just have to get the profile motif and compute the entropy for each column
        profile = profile_motifs(motifs)
        entropy_a = []
        global_entropy = 0.0
        (0..(profile.column_count-1)).each do |col_i|
            entropy = 0.0
            profile.column(col_i).each do |val|
                # puts val
                # puts (val*Math.log2(val))
                entropy += (val*Math.log2(val)) unless val == 0.0
                # puts entropy
            end
            entropy *= -1.0 unless entropy == 0.0
            entropy_a << entropy.round(4)
            global_entropy += entropy
        end    

        # return entropy_a
        return global_entropy.round(3)
    end # End of motifs_entropy function

    def greedy_motif_search(dna, k, t)
        # Our proposed greedy motif search algorithm, called GREEDYMOTIFSEARCH,
        # iteratively finds k-mers in the first string from Dna, second string from Dna, third string from Dna, etc.
        # After finding i - 1 k-mers Motifs in the first i - 1 strings of Dna, this algorithm constructs Profile(Motifs)
        # and selects the Profile-most probable k-mer from the i-th string based on this profile matrix.

        # GREEDYMOTIFSEARCH(Dna, k, t)
        #     BestMotifs ← motif matrix formed by first k-mers in each string
        #                   from Dna
        #     for each k-mer Motif in the first string from Dna
        #         Motif1 ← Motif
        #         for i = 2 to t
        #             form Profile from motifs Motif1, …, Motifi - 1
        #             Motifi ← Profile-most probable k-mer in the i-th string
        #                       in Dna
        #         Motifs ← (Motif1, …, Motift)
        #         if Score(Motifs) < Score(BestMotifs)
        #             BestMotifs ← Motifs
        #     output BestMotifs

        best_motifs = []
        dna.each do |dna_str|
            best_motifs << kmers(dna_str,k).first
        end
        kmers(dna[0], k).each do |kmer|
            motifs = [kmer]
            (1..(t-1)).each do |i|
                matrix = dna_matrix_from_kmers(motifs)
                profile = profile_motifs(matrix)
                motifs << profile_most_probable_kmer(dna[i], k, profile)
            end
            motifs_matrix = dna_matrix_from_kmers(motifs)
            best_motifs_matrix = dna_matrix_from_kmers(best_motifs)
            if score_motifs(motifs_matrix) < score_motifs(best_motifs_matrix)
                best_motifs = motifs.dup
            end
        end
        return best_motifs
    end # end of greedy_motif_search function

    def greedy_motif_search_laplace(dna, k, t)
        # Our proposed greedy motif search algorithm, called GREEDYMOTIFSEARCH,
        # iteratively finds k-mers in the first string from Dna, second string from Dna, third string from Dna, etc.
        # After finding i - 1 k-mers Motifs in the first i - 1 strings of Dna, this algorithm constructs Profile(Motifs)
        # and selects the Profile-most probable k-mer from the i-th string based on this profile matrix.

        # GREEDYMOTIFSEARCH(Dna, k, t)
        #     BestMotifs ← motif matrix formed by first k-mers in each string
        #                   from Dna
        #     for each k-mer Motif in the first string from Dna
        #         Motif1 ← Motif
        #         for i = 2 to t
        #             form Profile from motifs Motif1, …, Motifi - 1
        #             Motifi ← Profile-most probable k-mer in the i-th string
        #                       in Dna
        #         Motifs ← (Motif1, …, Motift)
        #         if Score(Motifs) < Score(BestMotifs)
        #             BestMotifs ← Motifs
        #     output BestMotifs

        best_motifs = []
        dna.each do |dna_str|
            best_motifs << kmers(dna_str,k).first
        end
        kmers(dna[0], k).each do |kmer|
            motifs = [kmer]
            (1..(t-1)).each do |i|
                matrix = dna_matrix_from_kmers(motifs)
                profile = profile_motifs_laplace(matrix)
                motifs << profile_most_probable_kmer(dna[i], k, profile)
            end
            motifs_matrix = dna_matrix_from_kmers(motifs)
            best_motifs_matrix = dna_matrix_from_kmers(best_motifs)
            if score_motifs(motifs_matrix) < score_motifs(best_motifs_matrix)
                best_motifs = motifs.dup
            end
        end
        return best_motifs
    end # end of greedy_motif_search_laplace function

    def gibbs_sampler(dna, k, t, n)
        # GIBBSSAMPLER(Dna, k, t, N)
        # randomly select k-mers Motifs = (Motif1, …, Motift) in each string
        #     from Dna
        # BestMotifs ← Motifs
        # for j ← 1 to N
        #     i ← Random(t)
        #     Profile ← profile matrix constructed from all strings in Motifs
        #                except for Motifi
        #     Motifi ← Profile-randomly generated k-mer in the i-th sequence
        #     if Score(Motifs) < Score(BestMotifs)
        #         BestMotifs ← Motifs
        # return BestMotifs

        # https://class.coursera.org/bioinformatics-002/forum/thread?thread_id=342#post-1398

        # Select a random kmer motif in each Dna string as your starter set of motifs (Maintain a record of which motif is picked for each Dna string)
        # Then, for each iteration, (run N iterations total)
            # Randomly pick one Dna string X to set aside (Let's say within X, the current motif is Y).
            # Formulate a new profile matrix (with pseudocounts), P, based on the motifs from the remaining Dna strings (i.e. excluding X).
            # Form a biased-die Bd, by computing the probability of each kmer motif in X based on P.
            # Roll Bd to find a new motif Z, to replace Y as your new selected motif for X.
            # Put X back into the set of Dna strings.

        motifs = []
        dna.each do |dna_str|
            dna_kmers = kmers(dna_str,k)
            rand_index = rand(dna_kmers.length)
            motifs << dna_kmers[rand_index]
        end
        best_motifs = motifs.dup    

        (1..n).each do |j|
            i = rand(t)
            # puts "Random number i: " + i.to_s
            motifs_except_i = []
            motifs_except_i = motifs.dup
            # puts "Motifs before deletion: " + motifs_except_i.to_s
            motifs_except_i.delete_at(i)
            # puts "Motifs after deletion: " + motifs_except_i.to_s
            motifs_except_i_matrix = dna_matrix_from_kmers(motifs_except_i)
            profile = profile_motifs_laplace(motifs_except_i_matrix)
            random_motif = profile_random_kmer_motif(dna[i], k, profile)
            # puts "random_motif: " + random_motif
            # puts motifs_except_i.join(" ")
            motifs_except_i.insert(i, random_motif)
            # puts motifs_except_i.join(" ")
            motifs = []
            motifs = motifs_except_i.dup
            motifs_matrix = dna_matrix_from_kmers(motifs)
            best_motifs_matrix = dna_matrix_from_kmers(best_motifs)
            if score_motifs(motifs_matrix) < score_motifs(best_motifs_matrix)
                best_motifs = motifs.dup
            end            
        end

        return best_motifs
    end # end of gibbs_sampler function

    def profile_random_kmer_motif(text, k, profile)
        # We have previously defined the notion of a Profile-most probable k-mer in a string.
        # We now define a Profile-randomly generated k-mer in a string Text. For each k-mer Pattern in Text,
        # compute the probability Pr(Pattern | Profile), resulting in n = |Text| - k + 1 probabilities (p1, …, pn).
        # These probabilities do not necessarily sum to 1, but we can still form the random number generator
        # Random(p1, …, pn) based on them. 

        # Profile-most rofile-randomly generated k-mer
        #     Input: A string Text, an integer k, and a 4 × k matrix Profile.
        #     Output: A Profile-random k-mer in Text.

        # puts profile
        # puts text
        # puts k
        kmers_arr = kmers(text,k)
        kmers_prob_arr = []
        kmers_arr.each do |kmer|
            kmers_prob_arr << dna_profile_probability(kmer, profile)
        end

        # prob_sum = 0.0
        # puts "prob_arr: " + kmers_prob_arr.join(" ")
        # kmers_prob_arr.each {|i| prob_sum += i}
        # puts "prob_sum: " + prob_sum.to_s
        prob_sum = kmers_prob_arr.inject(:+)
        # puts "inject sum: " + kmers_prob_arr.inject(:+).to_s
        normalized_probability_arr = kmers_prob_arr.map {|i| (i+0.0)/prob_sum}
        # puts "normalized_probability_arr: " + normalized_probability_arr.join(" ")
        # prob_sum = 0.0
        # normalized_probability_arr.each {|i| prob_sum += i}
        # puts "norm_prob_sum: " + prob_sum.to_s

        rand_num = biased_random_number(normalized_probability_arr)
        # puts kmers_arr
        # puts "Rand" + rand_num.to_s + "random"
        return kmers_arr[rand_num]
    end # End of profile_random_kmer_motif function

    def biased_random_number(prob_arr)
        # http://stackoverflow.com/questions/479236/how-do-i-simulate-biased-die-in-python
        rand_roll = rand()
        sum = 0 
        result = 0
        prob_arr.each do |i|
            sum += i
            if rand_roll < sum
                return result
            end
            result += 1
        end
        # puts "prob_arry:" + prob_arr.to_s
        # puts "rand_roll: " + rand_roll.to_s
        # puts prob_arr.join(" ")
        return result - 1
    end

    def randomized_motif_search_repeat(dna, k, t, n)
        best_motifs = randomized_motif_search(dna, k, t)
        (1..1000).each do |i|
            motifs = randomized_motif_search(dna, k, t)
            motifs_matrix = dna_matrix_from_kmers(motifs)
            best_motifs_matrix = dna_matrix_from_kmers(best_motifs)
            if score_motifs(motifs_matrix) < score_motifs(best_motifs_matrix)
                best_motifs = motifs.dup
            end
        end
        return best_motifs
    end # end of randomized_motif_search_repeat function

    def randomized_motif_search(dna, k, t)
        # RANDOMIZEDMOTIFSEARCH(Dna, k, t)
        #     randomly select k-mers Motifs = (Motif1, …, Motift) in each string
        #         from Dna
        #     BestMotifs ← Motifs
        #     while forever
        #         Profile ← Profile(Motifs)
        #         Motifs ← Motifs(Profile, Dna)
        #         if Score(Motifs) < Score(BestMotifs)
        #             BestMotifs ← Motifs
        #         else
        #             return BestMotifs        

        motifs = []
        dna.each do |dna_str|
            dna_kmers = kmers(dna_str,k)
            rand_index = rand(dna_kmers.length)
            motifs << dna_kmers[rand_index]
        end
        best_motifs = motifs.dup

        loop do
            matrix = dna_matrix_from_kmers(motifs)
            profile = profile_motifs_laplace(matrix)
            motifs = []
            motifs = motifs(profile, dna, k)
            motifs_matrix = dna_matrix_from_kmers(motifs)
            best_motifs_matrix = dna_matrix_from_kmers(best_motifs)
            if score_motifs(motifs_matrix) < score_motifs(best_motifs_matrix)
                best_motifs = motifs.dup
            else
                return best_motifs
            end
        end

    end # end of randomized_motif_search function

    def motifs(profile, dna, k)
        # We previously defined Profile(Motifs) as the profile matrix constructed from a collection of k-mers Motifs in Dna. 
        # Now, given a collection of strings Dna and an arbitrary 4 × k matrix Profile, we define Motifs(Profile,Dna)
        # as the collection of k-mers formed by the Profile-most probable k-mers in each sequence from Dna.
        # For example, consider the following Profile and Dna:
        motifs = []
        dna.each do |text|
            motifs << profile_most_probable_kmer(text.upcase, k, profile)
        end
        return motifs
    end

    def dna_matrix_from_kmers(kmers)
        # Given an array of DNA kmers, we have to construct a matrix of individual DNA elements
        # Ex:-
        # Input: ["AGC", "TCG"]
        # Output:
        # Matrix [["A","G","C"],
        #         ["T","C","G"]]

        rows = []
        kmers.each do |kmer|
            row = []
            kmer.chars.each do |dna_char|
                row << dna_char.upcase
            end
            rows << row
        end

        return Matrix.rows(rows)
    end # end of dna_matrix_from_kmers


    def profile_most_probable_kmer(text, k, profile)
        # Given a profile matrix Profile, we can evaluate the probability of every k-mer in a string Text
        # and find a Profile-most probable k-mer in Text, i.e., a k-mer that was most likely to have been
        # generated by Profile among all k-mers in Text. For example, ACGGGGATTACC is the Profile-most
        # probable 12-mer in GGTACGGGGATTACCT. Indeed, every other 12-mer in this string has probability 0.
        # In general, if there are multiple Profile-most probable k-mers in Text, then we select the first such k-mer occurring in Text.

        # Profile-most Probable k-mer Problem: Find a Profile-most probable k-mer in a string.
        #     Input: A string Text, an integer k, and a 4 × k matrix Profile.
        #     Output: A Profile-most probable k-mer in Text.

        # puts profile
        # puts text
        # puts k
        kmers_arr = kmers(text,k)
        mp_kmer = kmers_arr[0]
        # puts mp_kmer
        kmer_prob = dna_profile_probability(mp_kmer, profile)
        (1..(kmers_arr.length-1)).each do |i|
            kmer = kmers_arr[i]
            # puts kmer 
            prob = dna_profile_probability(kmer, profile)
            if prob > kmer_prob
                kmer_prob = prob
                mp_kmer = kmer
            end
        end

        return mp_kmer
    end # End of profile_most_probable_kmer function


    def dna_profile_probability(dna, profile)
        # Given the profile matrix, figure out the probability of a given DNA string
        # Profile
        #     A:  .2  .2   0   0   0   0  .9  .1  .1  .1  .3   0            
        #     C:  .1  .6   0   0   0   0   0  .4  .1  .2  .4  .6  
        #     G:   0   0   1   1  .9  .9  .1   0   0   0   0   0  
        #     T:  .7  .2   0   0  .1  .1   0  .5  .8  .7  .3  .4  

        # Consensus        
        # T   C   G   G   G   G   A   T   T   T   C   C              

        # Pr(TCGGGGATTTCC | Profile) = 0.7 · 0.6 · 1.0 · 1.0 · 0.9 · 0.9 · 0.9 · 0.5 · 0.8 · 0.7 · 0.4 · 0.6 = 0.0205753  

        # Lets have a row index map for quick access in the matrix
        dna_to_row = {"A" => 0, "C" => 1, "G" => 2, "T" => 3}

        prob = 1.0
        (0..(dna.length-1)).each do |j|
            # puts dna[j]
            # puts dna_to_row[dna[j]]
            prob *= profile.element(dna_to_row[dna[j]],j)
        end

        return prob
    end # end of dna_profile_probability function

    def profile_motifs(motifs)
        # motifs is a ruby Matrix
        # We can construct the 4 × k count matrix Count(Motifs) counting the number of occurrences of each nucleotide
        # in each column of the motif matrix; the (i, j)-th element of Count(Motifs) stores the number of times that
        # nucleotide i appears in column j of Motifs. We will further divide all of the elements in the count matrix by t,
        # the number of rows in Motifs. This results in a profile matrix P = Profile(Motifs) for which Pi,j is the frequency
        # of the i-th nucleotide in the j-th column of the motif matrix.
        # Note that the elements of any column of the profile matrix sum to 1

        # All we have to do is get the count matrix for the motifs and divide it up by the number of rows in the motifs Matrix
        count = count_motifs(motifs)
        
        # Adding 0.0 as a cheap trick to convert the division into a float
        return count.collect {|val| (val+0.0)/motifs.row_count}
    end

    def count_motifs(motifs)
        # motifs is a ruby Matrix
        # We can construct the 4 × k count matrix Count(Motifs) counting the number of occurrences of each nucleotide
        # in each column of the motif matrix; the (i, j)-th element of Count(Motifs) stores the number of times that
        # nucleotide i appears in column j of Motifs.

        count_arr = []
        (0..(motifs.column_count-1)).each do |col_i|
            col_array = motifs.column(col_i)
            # Each value of the count_col_array will correspond to the count of A, C, G, T respectively
            count_col_array = [0,0,0,0]
            col_array.each do |val|
                case val.upcase
                when "A"
                    count_col_array[0] += 1
                when "C"
                    count_col_array[1] += 1
                when "G"
                    count_col_array[2] += 1
                when "T"
                    count_col_array[3] += 1
                end
            end
            count_arr << count_col_array
        end

        # puts count_arr

        return Matrix.columns(count_arr)
    end # End of count_motifs function

   def profile_motifs_laplace(motifs)
        # motifs is a ruby Matrix
        # We can construct the 4 × k count matrix Count(Motifs) counting the number of occurrences of each nucleotide
        # in each column of the motif matrix; the (i, j)-th element of Count(Motifs) stores the number of times that
        # nucleotide i appears in column j of Motifs. We will further divide all of the elements in the count matrix by t,
        # the number of rows in Motifs. This results in a profile matrix P = Profile(Motifs) for which Pi,j is the frequency
        # of the i-th nucleotide in the j-th column of the motif matrix.
        # Note that the elements of any column of the profile matrix sum to 1

        # All we have to do is get the count matrix for the motifs and divide it up by the number of rows in the motifs Matrix
        count = count_motifs_laplace(motifs)
        
        # Adding 0.0 as a cheap trick to convert the division into a float
        return count.collect {|val| (val+0.0)/(motifs.row_count*2)}
    end # end of profile_motifs_laplace function

    def count_motifs_laplace(motifs)
        # motifs is a ruby Matrix
        # We can construct the 4 × k count matrix Count(Motifs) counting the number of occurrences of each nucleotide
        # in each column of the motif matrix; the (i, j)-th element of Count(Motifs) stores the number of times that
        # nucleotide i appears in column j of Motifs.

        count_arr = []
        (0..(motifs.column_count-1)).each do |col_i|
            col_array = motifs.column(col_i)
            # Each value of the count_col_array will correspond to the count of A, C, G, T respectively
            count_col_array = [1,1,1,1]
            col_array.each do |val|
                case val.upcase
                when "A"
                    count_col_array[0] += 1
                when "C"
                    count_col_array[1] += 1
                when "G"
                    count_col_array[2] += 1
                when "T"
                    count_col_array[3] += 1
                end
            end
            count_arr << count_col_array
        end

        # puts count_arr

        return Matrix.columns(count_arr)
    end # End of count_motifs_laplace function

    def score_motifs(motifs)
        # motifs is a ruby Matrix
        # To define scoring, consider t DNA sequences, each of length n, and select a k-mer 
        # from each sequence to form a collection Motifs, which we represent as a t × k motif matrix.
        # Our goal is to select k-mers resulting in the most “conserved” motif matrix, meaning the matrix with the 
        # most upper case letters (and thus the fewest number of lower case letters). Leaving aside the question of
        # how we select such k-mers, we will first focus on how to score the resulting motif matrices,
        # defining Score(Motifs) as the number of unpopular (lower case) letters in the motif matrix Motifs.
        # Our goal is to find a collection of k-mers that minimizes this score.

        # Motifs
        #    T   C   G   G   G   G   g   T   T   T   t   t           
        #    c   C   G   G   t   G   A   c   T   T   a   C
        #    a   C   G   G   G   G   A   T   T   T   t   C
        #    T   t   G   G   G   G   A   c   T   T   t   t
        #    a   a   G   G   G   G   A   c   T   T   C   C
        #    T   t   G   G   G   G   A   c   T   T   C   C
        #    T   C   G   G   G   G   A   T   T   c   a   t
        #    T   C   G   G   G   G   A   T   T   c   C   t
        #    T   a   G   G   G   G   A   a   c   T   a   C
        #    T   C   G   G   G   t   A   T   a   a   C   C

        # Score        
        #    3 + 4 + 0 + 0 + 1 + 1 + 1 + 5 + 2 + 3 + 6 + 4 = 30     

        # For now lets assume the input 'motifs' contain the upcase and lowercase encoded in it
        # If not we have to implement that in figuring out the score
        score = 0

        (0..(motifs.column_count-1)).each do |col_i|
            col_array = motifs.column(col_i)
            dna_freq = { "A" => 0, "C" => 0, "G" => 0, "T" => 0}
            col_array.each do |val|
                # score += 1 if val == val.downcase
                dna_freq[val.upcase] += 1
            end
            cnt_arr = dna_freq.values.sort
            count = 0
            cnt_arr.pop
            cnt_arr.each {|val| count += val}
            score += count
        end

        return score
    end # End of score_motifs function



    def all_else_dna(dna)
        # Input: A
        # Output: C G T
        # Input: C
        # Output: A G T
        dna_a = ["A", "C", "G", "T"]
        # puts dna_a.join(" ")
        dna_a.delete(dna)
        # puts dna_a.join(" ")
        return dna_a
    end # End of all_else_dna function

    def motif_enumeration(dna, k, d)
        # Given a collection of strings Dna and an integer d, a k-mer is a (k,d)-motif if it appears 
        # in every string from Dna with at most d mismatches. For example, the implanted 15-mer in the 
        # strings above represents a (15,4)-motif.

        # Implanted Motif Problem: Find all (k, d)-motifs in a collection of strings.
        #      Input: A collection of strings Dna, and integers k and d.
        #      Output: All (k, d)-motifs in Dna.

        # MOTIFENUMERATION(Dna, k, d)
        # Patterns ← an empty set
        # for each k-mer Pattern in Dna
        #     for each k-mer Pattern’ differing from Pattern by at most d
        #      mismatches
        #         if Pattern' appears in each string from Dna with at most d
        #      mismatches
        #             add Pattern' to Patterns
        # remove duplicates from Patterns
        # return Patterns        

        patterns = []
        (0..(dna[0].length-k)).each do |i|
            kmer = dna[0].slice(i,k)
            # puts "kmer:" + kmer
            iterative_neighbors(kmer, d).each do |pattern1|
                # puts "Pattern1: " + pattern1
                present_in_all = 1
                (1..(dna.length-1)).each do |j|
                    # puts "#{j}:" + approx_pattern_matching(pattern1, dna[j], d).join(" ")
                    # puts "#{j}:" + approx_pattern_matching(pattern1, dna[j], d).length.to_s
                    present_in_all = 0 if approx_pattern_matching(pattern1, dna[j], d).length == 0
                end
                if present_in_all == 1
                    patterns << pattern1 unless patterns.include?(pattern1)
                end
            end
        end

        return patterns
    end # End of motif_enumeration function

    def iterative_neighbors(pattern, d)
        # IterativeNeighbors(Pattern, d)
        # Neighborhood ← set consisting of single string Pattern
        # for j = 1 to d
        #     for each string Pattern’ in Neighborhood
        #         add ImmediateNeighbors(Pattern') to Neighborhood
        #         remove duplicates from Neighborhood
        # return Neighborhood

        neighborhood = [pattern]
        (1..d).each do |j|
            neighborhood.dup.each do |pattern1| 
                immediate_neighbors(pattern1).each do |neighbor|
                    neighborhood << neighbor
                end
            end
            neighborhood.uniq!
        end
        return neighborhood
    end # End of iterative_neighbors function


    def immediate_neighbors(pattern)
        # Input: ACG
        # Output: CCG TCG GCG AAG ATG AGG ACA ACC ACT ACG

        neighborhood = [pattern]
        (0..(pattern.length-1)).each do |i|
            symbol = pattern[i]
            # puts symbol
            # puts all_else_dna(symbol).join(" ")
            all_else_dna(symbol).each do |x|
                pattern_copy = pattern.dup
                pattern_copy[i] = x
                # puts pattern_copy
                neighborhood << pattern_copy
                # puts neighborhood.join(" ")
            end
        end
        return neighborhood
    end # End of immediate_neighbors function

    def peptide_mass(peptide)
        mass = 0
        amh = amino_acid_mass_hash()
        peptide.chars.each do |pep|
            mass += amh[pep]
        end
        return mass
    end # end of peptide_mass function

    def peptide_mass_with_amino_acids(peptide, amh)
        mass = 0
        peptide.chars.each do |pep|
            mass += amh[pep]
        end
        return mass
    end # end of peptide_mass_with_amino_acids function

    def expand_peptides(peptides)
        amh = amino_acid_mass_hash()
        exp_pep = []
        if peptides.empty?
            return amh.keys
        else 
            peptides.each do |pep|
                amh.keys.each do |aa|
                    exp_pep << (pep + aa)
                end
            end
        end 
        return exp_pep
    end # end of expand_peptides function

    def expand_peptides_with_amino_acids(peptides, amh)
        exp_pep = []
        if peptides.empty?
            return amh.keys
        else 
            peptides.each do |pep|
                amh.keys.each do |aa|
                    exp_pep << (pep + aa)
                end
            end
        end 
        return exp_pep
    end # end of expand_peptides_with_amino_acids function

    def peptide_consistency_with_spectrum?(peptide, spectrum)
        lin_spec_pep = linear_spectrum(peptide)

        lin_spec_pep.each do |pep|
            return false unless spectrum.include?(pep)
        end

        return true
    end # end of peptide_consistency_with_spectrum? function

    def peptide_to_mass(peptide)
        amh = amino_acid_mass_hash()
        mass_a = []
        peptide.chars.each {|p| mass_a << amh[p].to_s}
        mass_a.join("-")
    end # end of peptide_to_mass


    def peptide_to_mass_with_amino_acids(peptide, amh)
        mass_a = []
        peptide.chars.each {|p| mass_a << amh[p].to_s}
        mass_a.join("-")
    end # end of peptide_to_mass_with_amino_acids function

    def array_to_count_hash(arr)
        h = {}
        arr.each do |el|
            if h.has_key?(el)
                h[el] += 1
            else
                h[el] = 1
            end
        end
        return h
    end # end of array_to_count_hash function


    def linearpeptide_score_with_amino_acids(peptide, spectrum, amh)
        # Given a linear peptide Peptide and a spectrum Spectrum, we define LinearScore(Peptide, Spectrum) 
        # as the number of masses shared between Linearspectrum(Peptide) and Spectrum. 

        # Cyclopeptide Scoring Problem: Compute the score of a cyclic peptide against a spectrum.
        #      Input: An amino acid string Peptide and a collection of integers Spectrum. 
        #      Output: The score of Peptide against Spectrum, Score(Peptide, Spectrum).        

        theoretical_spectrum = linear_spectrum_with_amino_acids(peptide, amh)
        ts_h = array_to_count_hash(theoretical_spectrum)
        es_h = array_to_count_hash(spectrum)

        score = 0
        es_h.keys.each do |mass|
            if ts_h.has_key?(mass)
                # If the experimental spectrum has fewer instances of the same mass than theoritical spectrum use that one
                # to calculate the score but if experimental spectrum has more count of the same mass we should be using
                # the count in the theoritical spectrum
                if es_h[mass] <= ts_h[mass]
                    score += es_h[mass]
                else
                    score += ts_h[mass]
                end
            end
        end
        return score
    end # end of linearpeptide_score_with_amino_acids function

    def linearpeptide_score(peptide, spectrum)
        # Given a linear peptide Peptide and a spectrum Spectrum, we define LinearScore(Peptide, Spectrum) 
        # as the number of masses shared between Linearspectrum(Peptide) and Spectrum. 

        # Cyclopeptide Scoring Problem: Compute the score of a cyclic peptide against a spectrum.
        #      Input: An amino acid string Peptide and a collection of integers Spectrum. 
        #      Output: The score of Peptide against Spectrum, Score(Peptide, Spectrum).        

        theoretical_spectrum = linear_spectrum(peptide)
        ts_h = array_to_count_hash(theoretical_spectrum)
        es_h = array_to_count_hash(spectrum)

        score = 0
        es_h.keys.each do |mass|
            if ts_h.has_key?(mass)
                # If the experimental spectrum has fewer instances of the same mass than theoritical spectrum use that one
                # to calculate the score but if experimental spectrum has more count of the same mass we should be using
                # the count in the theoritical spectrum
                if es_h[mass] <= ts_h[mass]
                    score += es_h[mass]
                else
                    score += ts_h[mass]
                end
            end
        end
        return score
    end # end of linearpeptide_score function

    def cyclopeptide_score(peptide, spectrum)
        # Given a cyclic peptide Peptide and a spectrum Spectrum, we define Score(Peptide, Spectrum) 
        # as the number of masses shared between Cyclospectrum(Peptide) and Spectrum. 

        # Cyclopeptide Scoring Problem: Compute the score of a cyclic peptide against a spectrum.
        #      Input: An amino acid string Peptide and a collection of integers Spectrum. 
        #      Output: The score of Peptide against Spectrum, Score(Peptide, Spectrum).        

        theoretical_spectrum = cyclic_spectrum(peptide)
        ts_h = array_to_count_hash(theoretical_spectrum)
        es_h = array_to_count_hash(spectrum)

        score = 0
        es_h.keys.each do |mass|
            if ts_h.has_key?(mass)
                # If the experimental spectrum has fewer instances of the same mass than theoritical spectrum use that one
                # to calculate the score but if experimental spectrum has more count of the same mass we should be using
                # the count in the theoritical spectrum
                if es_h[mass] <= ts_h[mass]
                    score += es_h[mass]
                else
                    score += ts_h[mass]
                end
            end
        end
        return score
    end # end of cyclopeptide_score function

    def spectral_convolution(spectrum)
        # Convolution of a spectrum by taking all positive differences of masses in the spectrum.

        # Spectral Convolution Problem: Compute the convolution of a spectrum.
        #      Input: A collection of integers Spectrum.
        #      Output: The list of elements in the convolution of Spectrum. If an element has
        #      multiplicity k, it should appear exactly k times; you may return the elements in any order.
        spectrum.sort!
        convolution = []
        (1..(spectrum.length-1)).each do |i|
            (i-1).downto(0).each {|j| convolution << (spectrum[i] - spectrum[j])} 
        end
        convolution.delete(0)
        return convolution
    end # end of spectral_convolution function

    def leaderboard_trim(leader_board, spectrum, n)
        # Trim(Leaderboard, Spectrum, N, AminoAcid, AminoAcidMass)
        # for j ← 1 to |Leaderboard|
        #     Peptide ← j-th peptide in Leaderboard
        #     LinearScores(j) ← LinearScore(Peptide, Spectrum)
        # sort Leaderboard according to the decreasing order of scores in LinearScores
        # sort LinearScores in decreasing order
        # for j ← N + 1 to |Leaderboard|
        #     if LinearScores(j) < LinearScores(N)
        #         remove all peptides starting from the j-th peptide from Leaderboard
        #         return Leaderboard
        # return Leaderboard

        linear_scores_h = {}
        leader_board.each do |peptide|
            linear_scores_h[peptide] = linearpeptide_score(peptide, spectrum)
        end
        # puts linear_scores_h
        linear_scores_sorted_h = Hash[linear_scores_h.sort_by {|k,v| v}.reverse]
        # puts linear_scores_sorted_h
        (n..(linear_scores_sorted_h.length-1)).each do |j|
            if linear_scores_sorted_h.values[j] < linear_scores_sorted_h.values[n-1]
                leader_board = linear_scores_sorted_h.keys[0,j]
                return leader_board
            end
        end
        return leader_board
    end # end of leaderboard_trim function


    def leaderboard_trim_with_amino_acids(leader_board, spectrum, n, amh)
        # Trim(Leaderboard, Spectrum, N, AminoAcid, AminoAcidMass)
        # for j ← 1 to |Leaderboard|
        #     Peptide ← j-th peptide in Leaderboard
        #     LinearScores(j) ← LinearScore(Peptide, Spectrum)
        # sort Leaderboard according to the decreasing order of scores in LinearScores
        # sort LinearScores in decreasing order
        # for j ← N + 1 to |Leaderboard|
        #     if LinearScores(j) < LinearScores(N)
        #         remove all peptides starting from the j-th peptide from Leaderboard
        #         return Leaderboard
        # return Leaderboard

        linear_scores_h = {}
        leader_board.each do |peptide|
            linear_scores_h[peptide] = linearpeptide_score_with_amino_acids(peptide, spectrum, amh)
        end
        # puts linear_scores_h
        linear_scores_sorted_h = Hash[linear_scores_h.sort_by {|k,v| v}.reverse]
        # puts linear_scores_sorted_h
        (n..(linear_scores_sorted_h.length-1)).each do |j|
            if linear_scores_sorted_h.values[j] < linear_scores_sorted_h.values[n-1]
                leader_board = linear_scores_sorted_h.keys[0,j]
                return leader_board
            end
        end
        return leader_board
    end # end of leaderboard_trim_with_amino_acids function

    def frequent_elements_with_ties(spectrum, m)
        # Given an array of masses in the spectrum figure out the most 'm' frequent masses
        # including the ties

        freq_h = {}
        spectrum.each do |m|
            if freq_h[m]
                freq_h[m] += 1
            else
                freq_h[m] = 1
            end
        end

        freq_sorted_h = Hash[freq_h.sort_by{|k,v| v}.reverse]

        return_arr = freq_sorted_h.keys.slice(0, m)

        # puts freq_sorted_h

        (m..(freq_sorted_h.keys.length-1)).each do |i|
            if (freq_sorted_h.values[m-1] == freq_sorted_h.values[i])
                return_arr << freq_sorted_h.keys[i]
            elsif (freq_sorted_h.values[m-1] > freq_sorted_h.values[i])
                break
            end
        end
        # puts return_arr.to_s
        return return_arr
    end 

    def convolution_cyclopeptide_sequencing(m, n, spectrum)
        # We now have the outline for a new cyclopeptide sequencing algorithm. 
        # Given an experimental spectrum, we first compute the convolution of an experimental spectrum. 
        # We then select the M most frequent elements between 57 and 200 in the convolution to form an extended alphabet of amino acid masses. 
        # In order to be fair, we should include the top M elements of the convolution "with ties". 
        # Finally, we run the algorithm LeaderboardCyclopeptideSequencing, where the amino acid masses are restricted to this alphabet. 
        # We call this algorithm ConvolutionCyclopeptideSequencing.

        # CODE CHALLENGE: Implement ConvolutionCyclopeptideSequencing.
        #      Input: An integer M, an integer N, and a collection of (possibly repeated) integers Spectrum.
        #      Output: A cyclic peptide LeaderPeptide with amino acids taken only from the top M elements
        #      (and ties) of the convolution of Spectrum that fall between 57 and 200, and where the size
        #      of Leaderboard is restricted to the top N (and ties).
        
        # puts spectrum.join(" ")
        convol_spectrum = spectral_convolution(spectrum)
        convol_spectrum_57_200 = convol_spectrum.keep_if {|i| i >= 57 && i <= 200}
        # puts convol_spectrum.join(" ")
        aa_mass = frequent_elements_with_ties(convol_spectrum_57_200, m) 
        amh = amino_acid_mass_mapping(aa_mass)    
        puts amh
        return leaderboard_cyclopeptide_sequencing_with_amino_acids(spectrum, n , amh)
    end # end of convolution_cyclopeptide_sequencing

    def amino_acid_mass_mapping(aa_mass)
        amh = {}
        aa_mass.each do |m|
            amh[m.chr] = m
        end
        return amh
    end


    def leaderboard_cyclopeptide_sequencing_with_amino_acids(spectrum, n, amh)
        # To be fair, a cut should include anyone who is tied with the Nth-place competitor.
        # Thus, Leaderboard should be trimmed down to the “N highest-scoring linear peptides including ties”, 
        # which may include more than N peptides. Given a list of peptides Leaderboard, a spectrum Spectrum, 
        # and an integer N, define Trim(Leaderboard, Spectrum, N) as the collection of the top N highest-scoring 
        # linear peptides in Leaderboard (including ties) with respect to Spectrum.

        # Note that Score(Peptide, Spectrum) currently only scores Peptide against Spectrum if Peptide is cyclic. 
        # However, to generalize this scoring function when Peptide is linear, we simply exclude those subpeptides
        # of Peptide that wrap around the end of the string, resulting in a function LinearScore(Peptide, Spectrum). 
        # For example, if Spectrum is the experimental spectrum of NQEL, then you can verify that LinearScore(Peptide, Spectrum) = 8.

        # LEADERBOARDCYCLOPEPTIDESEQUENCING(Spectrum, N)
        #     Leaderboard ← {empty peptide}
        #     LeaderPeptide ← empty peptide
        #     while Leaderboard is non-empty
        #         Leaderboard ← Expand(Leaderboard)
        #         for each Peptide in Leaderboard
        #             if Mass(Peptide) = ParentMass(Spectrum)
        #                 if Score(Peptide, Spectrum) > Score(LeaderPeptide, Spectrum)
        #                     LeaderPeptide ← Peptide
        #             else if Mass(Peptide) > ParentMass(Spectrum)
        #                 remove Peptide from Leaderboard
        #         Leaderboard ← Trim(Leaderboard, Spectrum, N)
        #     output LeaderPeptide      

        leader_board = []
        leader_peptide = ""
        parent_mass = spectrum[-1]
        loop do
            leader_board = expand_peptides_with_amino_acids(leader_board, amh)
            # puts leader_board.join(" ")
            leader_board.dup.each do |peptide|
                # puts parent_mass.to_s + " " + peptide_mass(peptide).to_s
                if peptide_mass_with_amino_acids(peptide, amh)  == parent_mass
                    if linearpeptide_score_with_amino_acids(peptide, spectrum, amh) > linearpeptide_score_with_amino_acids(leader_peptide, spectrum, amh)
                        leader_peptide = peptide
                    end
                elsif peptide_mass_with_amino_acids(peptide, amh) > parent_mass
                    leader_board.delete(peptide)
                end
            end
            # puts leader_board.join(" ")
            leader_board = leaderboard_trim_with_amino_acids(leader_board, spectrum, n, amh) 
            break if leader_board.empty?
        end
        return peptide_to_mass_with_amino_acids(leader_peptide, amh)
    end # end of leaderboard_cyclopeptide_sequencing_with_amino_acids function


    def leaderboard_cyclopeptide_sequencing(spectrum, n)
        # To be fair, a cut should include anyone who is tied with the Nth-place competitor.
        # Thus, Leaderboard should be trimmed down to the “N highest-scoring linear peptides including ties”, 
        # which may include more than N peptides. Given a list of peptides Leaderboard, a spectrum Spectrum, 
        # and an integer N, define Trim(Leaderboard, Spectrum, N) as the collection of the top N highest-scoring 
        # linear peptides in Leaderboard (including ties) with respect to Spectrum.

        # Note that Score(Peptide, Spectrum) currently only scores Peptide against Spectrum if Peptide is cyclic. 
        # However, to generalize this scoring function when Peptide is linear, we simply exclude those subpeptides
        # of Peptide that wrap around the end of the string, resulting in a function LinearScore(Peptide, Spectrum). 
        # For example, if Spectrum is the experimental spectrum of NQEL, then you can verify that LinearScore(Peptide, Spectrum) = 8.

        # LEADERBOARDCYCLOPEPTIDESEQUENCING(Spectrum, N)
        #     Leaderboard ← {empty peptide}
        #     LeaderPeptide ← empty peptide
        #     while Leaderboard is non-empty
        #         Leaderboard ← Expand(Leaderboard)
        #         for each Peptide in Leaderboard
        #             if Mass(Peptide) = ParentMass(Spectrum)
        #                 if Score(Peptide, Spectrum) > Score(LeaderPeptide, Spectrum)
        #                     LeaderPeptide ← Peptide
        #             else if Mass(Peptide) > ParentMass(Spectrum)
        #                 remove Peptide from Leaderboard
        #         Leaderboard ← Trim(Leaderboard, Spectrum, N)
        #     output LeaderPeptide      

        leader_board = []
        leader_peptide = ""
        parent_mass = spectrum[-1]
        loop do
            leader_board = expand_peptides(leader_board)
            # puts leader_board.join(" ")
            leader_board.dup.each do |peptide|
                # puts parent_mass.to_s + " " + peptide_mass(peptide).to_s
                if peptide_mass(peptide)  == parent_mass
                    if linearpeptide_score(peptide, spectrum) > linearpeptide_score(leader_peptide, spectrum)
                        leader_peptide = peptide
                    end
                elsif peptide_mass(peptide) > parent_mass
                    leader_board.delete(peptide)
                end
            end
            # puts leader_board.join(" ")
            leader_board = leaderboard_trim(leader_board, spectrum, n) 
            break if leader_board.empty?
        end
        return peptide_to_mass(leader_peptide)
    end # end of leaderboard_cyclopeptide_sequencing function

    def cyclopeptide_sequencing(spectrum)
        # CYCLOPEPTIDESEQUENCING(Spectrum)
        #     Peptides ← {empty peptide}
        #     while Peptides is nonempty
        #         Peptides ← Expand(Peptides)
        #         for each peptide Peptide in Peptides
        #             if Mass(Peptide) = ParentMass(Spectrum)
        #                 if Cyclospectrum(Peptide) = Spectrum
        #                     output Peptide
        #                 remove Peptide from Peptides
        #             else if Peptide is not consistent with Spectrum
        #                 remove Peptide from Peptides

        # What about the branching step? Given the current collection of linear peptides Peptides, 
        # define Expand(Peptides) as a new collection containing all possible extensions of peptides in 
        # Peptides by a single amino acid mass
        peptides = []
        parent_mass = spectrum[-1]
        output_peptides = []
        loop do
            expanded_peptides = expand_peptides(peptides)
            # puts peptides
            delete_peptides = []
            expanded_peptides.each do |peptide|
                # puts peptide
                if peptide_mass(peptide) == parent_mass
                    cyclic_spec_pep = cyclic_spectrum(peptide)
                    # puts peptide
                    # puts cyclic_spec_pep.join(" ")
                    if cyclic_spec_pep == spectrum        
                        # puts "output: " + peptide_to_mass(peptide)
                        output_peptides << peptide_to_mass(peptide) 
                    end
                    delete_peptides << peptide
                else
                    unless peptide_consistency_with_spectrum?(peptide, spectrum)
                        delete_peptides << peptide
                    end
                end
            end
            delete_peptides.uniq.each {|p| expanded_peptides.delete(p) }

            break if expanded_peptides.empty?
            peptides = []
            expanded_peptides.each {|p| peptides << p}
        end

        output_peptides.uniq
    end # end of cyclopeptide_sequencing function


    def cyclic_spectrum(peptide)
        # If Peptide represents a cyclic peptide instead, then the masses in its theoretical spectrum 
        # can be divided into those found by LinearSpectrum and those corresponding to subpeptides 
        # wrapping around the end of Peptide. Furthermore, each such subpeptide has mass equal to the
        # difference between Mass(Peptide) and a subpeptide mass identified by LinearSpectrum. 
        # For example, when Peptide = NQEL, the mass of LN is equal to Mass(NQEL) minus the mass of QE, 
        # or 484 − 257 = 227. Thus, we can generate a cyclic spectrum by making only a small modification to the pseudocode of LinearSpectrum.

        # CyclicSpectrum(Peptide, AminoAcid, AminoAcidMass)
        #     PrefixMass(0) ← 0
        #     for i ← 1 to |Peptide|
        #         for j ← 1 to 20
        #             if AminoAcid(j) =  i-th amino acid in Peptide
        #                 PrefixMass(i) ← PrefixMass(i − 1) + AminoAcidMass(j)
        #     peptideMass ← PrefixMass(|Peptide|)
        #     CyclicSpectrum ← a list consisting of the single integer 0
        #     for i ← 0 to |Peptide| − 1
        #         for j ← i + 1 to |Peptide|
        #             add PrefixMass(j) − PrefixMass(i) to CyclicSpectrum
        #             if i > 0 and j < |Peptide|
        #                 add peptideMass - (PrefixMass(j) − PrefixMass(i)) to CyclicSpectrum
        #     return sorted list CyclicSpectrum

        amh = amino_acid_mass_hash()
        prefix_mass = [0]
        (1..peptide.length).each do |i|
            prefix_mass[i] = prefix_mass[i-1] + amh[peptide[i-1]]
        end

        peptide_mass = prefix_mass[peptide.length]
        spectrum = [0]
        (0..(peptide.length-1)).each do |i|
            ((i+1)..peptide.length).each do |j|
                spectrum << (prefix_mass[j] - prefix_mass[i])
                if i > 0 && j < peptide.length
                    spectrum << (peptide_mass - (prefix_mass[j] - prefix_mass[i]))
                end
            end
        end

        spectrum.sort
    end # end of cyclic_spectrum function

    def subpeptides_in_linear_peptide(n)
        # How many subpeptides does a linear peptide of given length n have? (Include the empty peptide and the entire peptide.)
        #      Input: An integer n.
        #      Output: The number of subpeptides of a linear peptide of length n.

        return ( ((n*(n+1))/2)  + 1 )
    end # end of subpeptides_in_linear_peptide

    def linear_spectrum_with_amino_acids(peptide, amh)
        # LinearSpectrum(Peptide, AminoAcid, AminoAcidMass)
        #     PrefixMass(0) ← 0
        #     for i ← 1 to |Peptide|
        #         for j ← 1 to 20
        #             if AminoAcid(j) =  i-th amino acid in Peptide
        #                 PrefixMass(i) ← PrefixMass(i − 1) + AminoAcidMass(j)
        #     LinearSpectrum ← a list consisting of the single integer 0
        #     for i ← 0 to |Peptide| − 1
        #         for j ← i + 1 to |Peptide|
        #             add PrefixMass(j) − PrefixMass(i) to LinearSpectrum
        #     return sorted list LinearSpectrum

        # CODE CHALLENGE: Implement LinearSpectrum.
        #      Input: An amino acid string Peptide.
        #      Output: The linear spectrum of Peptide.

        prefix_mass = [0]
        (1..peptide.length).each do |i|
            prefix_mass[i] = prefix_mass[i-1] + amh[peptide[i-1]]
        end

        spectrum = [0]
        (0..(peptide.length-1)).each do |i|
            ((i+1)..peptide.length).each do |j|
                spectrum << (prefix_mass[j] - prefix_mass[i])
            end
        end

        spectrum.sort
    end # end of linear_spectrum_with_amino_acids function


    def linear_spectrum(peptide)
        # LinearSpectrum(Peptide, AminoAcid, AminoAcidMass)
        #     PrefixMass(0) ← 0
        #     for i ← 1 to |Peptide|
        #         for j ← 1 to 20
        #             if AminoAcid(j) =  i-th amino acid in Peptide
        #                 PrefixMass(i) ← PrefixMass(i − 1) + AminoAcidMass(j)
        #     LinearSpectrum ← a list consisting of the single integer 0
        #     for i ← 0 to |Peptide| − 1
        #         for j ← i + 1 to |Peptide|
        #             add PrefixMass(j) − PrefixMass(i) to LinearSpectrum
        #     return sorted list LinearSpectrum

        # CODE CHALLENGE: Implement LinearSpectrum.
        #      Input: An amino acid string Peptide.
        #      Output: The linear spectrum of Peptide.

        amh = amino_acid_mass_hash()
        prefix_mass = [0]
        (1..peptide.length).each do |i|
            prefix_mass[i] = prefix_mass[i-1] + amh[peptide[i-1]]
        end

        spectrum = [0]
        (0..(peptide.length-1)).each do |i|
            ((i+1)..peptide.length).each do |j|
                spectrum << (prefix_mass[j] - prefix_mass[i])
            end
        end

        spectrum.sort
    end # end of linear_spectrum function

    def rna_peptide_match?(rna, peptide)
        codon_length = 3
        i = 0
        j = 0
        while (codon = rna.slice(i, codon_length)) do
            # puts rna_to_codon_hash[codon.to_sym] + " " + peptide[j]
            return false if rna_to_codon_hash[codon.to_sym] != peptide[j]
            i += codon_length
            j += 1
        end
        # puts "true"
        return true
    end # end of rna_peptide_match

    def reverse_complement(txt)
        comp = []
        txt.chars.each do |ch|
            comp << "A" if ch == 'T'
            comp << "T" if ch == 'A'
            comp << "G" if ch == 'C'
            comp << "C" if ch == 'G'
        end
        comp.join.reverse
    end # end of reverse_complement function

    def encode_peptide(genome, amino_acid)
        # We say that a DNA string Pattern encodes an amino acid string Peptide if the 
        # RNA string transcribed from either Pattern or its reverse complement Pattern translates into Peptide. 
        # For example, the DNA string GAAACT is transcribed into GAAACU and translated into ET. 
        # The reverse complement of this DNA string, AGTTTC, is transcribed into AGUUUC and translated into SF. 
        # Thus, GAAACT encodes both ET and SF.”

        # Peptide Encoding Problem: Find substrings of a genome encoding a given amino acid sequence.
        #      Input: A DNA string Text, an amino acid string Peptide, and the array GeneticCode.
        #      Output: All substrings of Text encoding Peptide (if any such substrings exist).        

        match_dna = []
        codon_length = 3
        rna_seq_length = (amino_acid.length * codon_length)
        
        rna = dna_to_rna(genome)
        i = 0
        while (rna_seq = rna.slice(i, rna_seq_length )) do
            # puts rna_seq
            match_dna <<  rna_to_dna(rna_seq) if rna_peptide_match?(rna_seq, amino_acid)
            i += 1
        end

        rna = dna_to_rna(reverse_complement(genome))
        # puts rna
        i = 0
        while (rna_seq = rna.slice(i, rna_seq_length )) do
            # puts rna_seq
            match_dna <<  reverse_complement(rna_to_dna(rna_seq)) if rna_peptide_match?(rna_seq, amino_acid)
            i += 1
        end

        match_dna
    end # end of encode_peptide function

    def rna_to_amino_acid(rna)
        # Protein Translation Problem: Translate an RNA string into an amino acid string.
        #      Input: An RNA string Pattern and the array GeneticCode.
        #      Output: The translation of Pattern into an amino acid string Peptide.        

        r_to_c_h = rna_to_codon_hash
        # puts r_to_c_h
        i = 0
        codon_length = 3
        amino_acid = ""
        while (codon = rna.slice(i, codon_length)) do
            # puts codon
            # puts r_to_c_h[codon.to_sym]
            break if codon.empty?
            amino_acid += r_to_c_h[codon.to_sym].to_s
            i += codon_length
        end
        return amino_acid
    end # end of rna_to_amino_acid function

    def d_neighbors_of_k_mer(pattern, k, d)
        # The d-neighborhood of the k-mer Pattern is the collection of all k-mers 
        # that are at most Hamming distance d from Pattern.
        #  How many 4-mers are in the 3-neighborhood of Pattern = TAGC? 
        # Note that the d-neighborhood of Pattern includes Pattern.

        genome_alphabet = ["A","G","C","T"]

        count = 0
        genome_alphabet.repeated_permutation(k).to_a.each do |kmer|
            count += 1 if hamming_distance(kmer.join(""), pattern) <= d
        end
        return count
    end # end of d_neighbors_of_k_mer

    def approx_pattern_count(text, pattern, d)
        # ApproximatePatternCount(Text, Pattern, d)
        #     count ← 0
        #     for i ← 0 to |Text| − |Pattern|
        #         Pattern′ ← Text(i , |Pattern|)
        #         if HammingDistance(Pattern, Pattern′) ≤ d
        #             count ← count + 1
        #     return count        

        # We are doing pretty much the same thing in approx pattern matching, we can just
        # get the array of indexes of matched patterns and count that

        return approx_pattern_matching(pattern, text, d).length
    end # end of approx_pattern_count


    def approx_pattern_matching(pattern, text, d)
        # We say that a k-mer Pattern appears as a substring of Text with at most d mismatches if there is 
        # some k-mer substring Pattern' of Text having d or fewer mismatches with Pattern, i.e., 
        # HammingDistance(Pattern, Pattern') ≤ d. Our observation that a DnaA box may appear with slight variations 
        # leads to the following generalization of the Pattern Matching Problem.

        # Approximate Pattern Matching Problem: Find all approximate occurrences of a pattern in a string.
        #      Input: Strings Pattern and Text along with an integer d.
        #      Output: All starting positions where Pattern appears as a substring of Text with at most d mismatches.

        # Exact pattern matching algo
        # match_indexes = []
        # search_start_pos = 0
        # while index = genome.index(pattern, search_start_pos)
        #     match_indexes << index
        #     # puts index
        #     search_start_pos = index + 1
        # end
        # return match_indexes        

        match_indexes = []
        search_start_pos = 0
        (0..((text.length - pattern.length))).each do |i|
            match_indexes << i if (hamming_distance(text.slice(i, pattern.length), pattern) <= d)
        end
        return match_indexes
    end # end of approx_pattern_matching function


    def hamming_distance(p, q)
        # Hamming Distance Problem: Compute the Hamming distance between two strings.
        #      Input: Two strings of equal length.
        #      Output: The Hamming distance between these strings.

        # We say that position i in k-mers p1 … pk and q1 … qk is a mismatch if p(i) ≠ q(i). 
        # For example, CGAAT and CGGAC have two mismatches. 
        # The number of mismatches between strings p and q is called the Hamming distance 
        # between these strings and is denoted HammingDistance(p, q).
        distance = 0
        (0..(p.length - 1)).each {|i| distance += 1 if p[i] != q[i] }                         

        return distance
    end # end of hamming_distance function


    def minimum_skew(genome)
        # Minimum Skew Problem: Find a position in a genome minimizing the skew.
        #      Input: A DNA string Genome.
        #      Output: All integer(s) i minimizing Skewi (Genome) among all values of i (from 0 to |Genome|).        

        skew_vals = skew(genome)
        min_skew = skew_vals[0]
        min_skew_indices = []
        min_skew_indices << 0
        (1..(skew_vals.length-1)).each do |i|
            if skew_vals[i] < min_skew
                min_skew = skew_vals[i]
                min_skew_indices = [i]
            elsif skew_vals[i] == min_skew
                min_skew_indices << i
            end
        end

        return min_skew_indices
    end # end of minimum_skew function


    def skew(genome)
        # Note that we can compute Skewi+1(Genome) from Skewi(Genome) according to the nucleotide in position i of Genome. 
        # If this nucleotide is G, then Skewi+1(Genome) = Skewi(Genome) + 1; 
        # if this nucleotide is C, then Skewi+1(Genome)= Skewi(Genome) – 1; otherwise, Skewi+1(Genome) = Skewi(Genome). 

        skew = []
        g_minus_c = 0
        skew << g_minus_c
        genome.chars.each do |dna|
            if dna == "C"
                g_minus_c = g_minus_c - 1 
            elsif dna == "G"
                g_minus_c = g_minus_c + 1 
            end
            skew << g_minus_c
        end       

        return skew
    end # end of skew function

    def number_to_pattern(index, k)
        # NumberToPattern(index, k)
        #     if k = 1
        #         return NumberToSymbol(index)
        #     prefixIndex ← Quotient(index, 4)
        #     r ← Remainder(index, 4)
        #     PrefixPattern ← NumberToPattern(prefixIndex, k − 1)
        #     symbol ← NumberToSymbol(r)
        #     return concatenation of PrefixPattern with symbol        


        # Input: Integers index and k.
        # Output: The string NumberToPattern(index, k). 
        # Sample Input:
        # 45
        # 4
        # Sample Output:
        # AGTC        

        return number_to_symbol(index) if k == 1
        prefix_index = index/4
        r = index%4
        prefix_pattern = number_to_pattern(prefix_index, k - 1)
        symbol = number_to_symbol(r)
        return prefix_pattern + symbol
    end # end of number_to_pattern


    def number_to_symbol(number)
        if number.to_i == 0
            return "A"
        elsif number.to_i == 1
            return "C"
        elsif number.to_i == 2
            return "G"
        elsif number.to_i == 3
            return "T"
        end
                
        return ""
    end # end of number_to_symbol

    def symbol_to_number(symbol)
        if symbol == "A"
            return 0
        elsif symbol == "C"
            return 1
        elsif symbol == "G"
            return 2
        elsif symbol == "T"
            return 3
        end

        return -1        
    end # end of symbol_tonumber function


    def pattern_to_number(pattern)
       # PatternToNumber(Pattern)
       #      if Pattern contains no symbols
       #          return 0
       #      symbol ← LastSymbol(Pattern)
       #      remove LastSymbol(Pattern) from Pattern
       #      return 4 · PatternToNumber(Pattern) + SymbolToNumber(symbol)

       # CODE CHALLENGE: Implement PatternToNumber.
       # Input: A DNA string Pattern.
       # Output: The integer PatternToNumber(Pattern).

        # Sample Input:
        # AGT
        # Sample Output:
        # 11       

        return 0 if pattern.empty?
        symbol = pattern.slice(pattern.length - 1)
        pattern.chop!
        return (4 * pattern_to_number(pattern)) + symbol_to_number(symbol)

    end # end of pattern_to_number function


    def faster_frequent_words(text, k)
       # FasterFrequentWords(Text , k)
       #      FrequentPatterns ← an empty set
       #      FrequencyArray ← ComputingFrequencies(Text, k)
       #      maxCount ← maximal value in FrequencyArray
       #      for i ←0 to 4**k − 1
       #          if FrequencyArray(i) = maxCount
       #              Pattern ← NumberToPattern(i, k)
       #              add Pattern to the set FrequentPatterns
       #      return FrequentPatterns

       #  Input: A string Text and an integer k.
       #  Output: All most frequent k-mers in Text.

       # Sample Input:
       #      ACGTTGCATGTCGCATGATGCATGAGAGCT
       #      4

       # Sample Output:
       #      CATG GCAT

       frequent_patterns = []
       frequency_array = frequency_array(text, k)
       max_count = 0
       frequency_array.each {|count| max_count = count if count > max_count }
       (0..(4**k - 1)).each do |i|
            if frequency_array[i] == max_count
                pattern = number_to_pattern(i,k)
                frequent_patterns << pattern
            end
       end

       return frequent_patterns
    end # end of faster_frequent_words


    def frequency_array(text, k)
        # ComputingFrequencies(Text , k)
        #     for i ← 0 to 4**k − 1
        #         FrequencyArray(i) ← 0
        #     for i ← 0 to |Text| − k
        #         Pattern ← Text(i, k)
        #         j ← PatternToNumber(Pattern)
        #         FrequencyArray(j) ← FrequencyArray(j) + 1
        #     return FrequencyArray        

        # CODE CHALLENGE: Implement ComputingFrequencies to generate a frequency array.
        #     Input: A DNA string Text followed by an integer k.
        #     Output: FrequencyArray(Text).        

        # Sample Input:
        # ACGCGGCTCTGAAA
        # 2
        # Sample Output:
        # 2 1 0 0 0 0 2 2 1 2 1 0 0 1 1 0        
        frequency_array = []
        (0..(4**k - 1)).each {|i| frequency_array[i]  = 0 }

        (0..(text.length - k)).each do |i|
            pattern = text.slice(i,k)
            j = pattern_to_number(pattern)
            frequency_array[j] = frequency_array[j] +  1
        end

        return frequency_array
    end

    def find_clump(genome, k, l, t)
        # Clump Finding Problem: Find patterns forming clumps in a string.
        #      Input: A string Genome, and integers k, L, and t.
        #      Output: All distinct k-mers forming (L, t)-clumps in Genome.

        # Sample Input:
        #      CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA
        #      5 50 4

        # Sample Output:
        #      CGACA GAAGA        

        # ClumpFinding(Genome, k, t, L)
        #     FrequentPatterns ← an empty set
        #     for i ←0 to 4k − 1
        #         Clump(i) ← 0
        #     for i ← 0 to |Genome| − L
        #         Text ← the string of length L starting at position i in Genome 
        #         FrequencyArray ← ComputingFrequencies(Text, k)
        #         for j ← 0 to 4k − 1
        #             if FrequencyArray(j) ≥ t
        #                 Clump(j) ← 1
        #     for i ← 0 to 4k − 1
        #         if Clump(i) = 1
        #             Pattern ← NumberToPattern(i, k)
        #             add Pattern to the set FrequentPatterns
        #     return FrequentPatterns


        # BetterClumpFinding(Genome, k, t, L)
        #     FrequentPatterns ← an empty set
        #     for i ←0 to 4**k − 1
        #         Clump(i) ← 0
        #     Text ← Genome(0, L)
        #     FrequencyArray ← ComputingFrequencies(Text, k)
        #     for i ← 0 to 4**k − 1
        #         if FrequencyArray(i) ≥ t
        #             Clump(i) ← 1
        #     for i ← 1 to |Genome| − L
        #         FirstPattern ← Genome(i − 1, k)
        #         j ← PatternToNumber(FirstPattern)
        #         FrequencyArray(j) ← FrequencyArray(j) − 1
        #         LastPattern ← Genome(i + L − k, k)
        #         j ← PatternToNumber(LastPattern)
        #         FrequencyArray(j) ← FrequencyArray(j) + 1
        #         if FrequencyArray(j) ≥ t
        #             Clump(j) ← 1
        #     for i ← 0 to 4**k − 1
        #         if Clump(i) = 1
        #             Pattern ← NumberToPattern(i, k)
        #             add Pattern to the set FrequentPatterns
        #     return FrequentPatterns

        frequent_patterns = []
        clump = []
        (0..(4**k -1)).each {|i| clump[i] = 0}
        text = genome.slice(0,l)
        frequency_array = frequency_array(text, k)
        (0..(4**k - 1)).each {|i| clump[i] = 1 if frequency_array[i] >= t}
        (1..(genome.length - l)).each do |i|
          first_pattern = genome.slice(i-1, k)
          j = pattern_to_number(first_pattern)
          frequency_array[j] = frequency_array[j] - 1
          last_pattern = genome.slice(i+l-k, k)
          j = pattern_to_number(last_pattern)
          frequency_array[j] = frequency_array[j] + 1
          clump[j] = 1 if frequency_array[j] >= t
        end
        (0..(4**k-1)).each do |i|
          if clump[i] == 1
            pattern = number_to_pattern(i, k)
            frequent_patterns << pattern
          end
        end
        
        return frequent_patterns
    end # end of find_clump function



    def pattern_matching(pattern, genome)
        # Pattern Matching Problem: Find all occurrences of a pattern in a string.
        #      Input: Two strings, Pattern and Genome.
        #      Output: All starting positions where Pattern appears as a substring of Genome.

        # Sample Input:
        #      ATAT
        #      GATATATGCATATACTT

        # Sample Output:
        #      1 3 9

        match_indexes = []
        search_start_pos = 0
        while index = genome.index(pattern, search_start_pos)
            match_indexes << index
            # puts index
            search_start_pos = index + 1
        end
        return match_indexes
    end # end of pattern_matching function

end # end of class BioInfoAlgos

# 