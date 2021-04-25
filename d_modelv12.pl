# Version 6 of Dispersal Model:
# Migration Rate Variables are in place for high, mid and low dispersal phenotypes.  
# Figured out Habitat 7 (sort of)
# Add output file that has allele frequency information


############################################################
######################## Parameters ########################
############################################################
use constant RANDBITS => 15; 
use constant RAND_MAX => 2**RANDBITS;

srand(time ^ ($$ + ($$ << 15)) ); # Seeds the PRNG

#$| = 1;
#select((select(FH9), $|=1)[0]);
#open (FH9, ">>test.csv");
#select((select(FH9), $|=1)[0]);
#print(FH9 "test");
### Logistics of model 
$runs = 1000; # Number of runs
$invgen = 0; # Generation when individuals can start migrating
$endgen = 10000; 


### Model Parameters 
$carrying_capacity = 500; # Carrying capacity 
$mutrt = 0.000002; # Nuclear mutation rate   
#$mutrt = 0.02; # Nuclear mutation rate  
#$cancost_rand = 6.05; #Used to seed the random number generator to determine the $cancost value each time the mother has a cannibal phenotype --cancost = 1 below 
#$dead = 0.25; #the probability of expressing the cannibal phenotype if an individual has the cannibal genotype
$num_cannibal_offspring = 1; #the number of additional offspring that a cannibal female can produce each mating 

$chance_cannibalism_cost = 0;  #probability that a cannibal mother will eat others of her population 
$chance_cannibalism_offspring = 0;  #probability that a cannibal mother will have additional offspring 


# Species A and B differ in migration rates 
$migrtH = 0.0005;  # Dispersal rate phenotype A
$migrtL = 0.00005;  # Dispersal rate phenotype B
#$migrtH = 0.00;  # Dispersal rate phenotype A
#$migrtL = 0.000;  # Dispersal rate phenotype B
$migrtM = ($migrtH + $migrtL) / 2;  # Intermediate Dispersal Rate for heterozygoats
# open(FH2, ">posterT.csv"); 

### Reproductive Control 
# Used to determine the number of offspring each female will have (average offspring per female is 2.05 1 (2.2) )

$Rexp = 20; 
$master_R = (2.05); 

$femrep = 1; # Control of reproductive output; 1=females, 0=both sexes
$offSR = 0.5; # Average proportion of females in offspring 


### Fitness Costs 
# CURRENTLY NOT IN USE
$flatFC = 0.0; # Flat fitness cost (probability of offspring death when mother is a non-cannibal) 
$canFC0 = $canFC1 = $canFC2 = 0.0; # Cannibal fitness cost (probability of offspring death when cannibal mother is of non-cannibal genotype or cannibal genotype) 


@all_loop_costs = (0.15); # cost to the population from having a cannibalistic individual
@all_loop_benefits = (0.075); # benefit to individual for being a cannibal
#@all_loop_costs = (0.0); # cost to the population from having a cannibalistic individual
#@all_loop_benefits = (0.2); # benefit to individual for being a cannibal




# Immigration Rates--don't mess with 
$staySameDisH = 1 - $migrtH;
$moveNewDisH = (1 - $staySameDisH) / 9;

$staySameDisM = 1 - $migrtM;
$moveNewDisM = (1 - $staySameDisM) / 9;  

$staySameDisL = 1 - $migrtL;
$moveNewDisL = (1 - $staySameDisL) / 9;  


# File Name
$filename = "$migrtH$migrtL";
$sumname = "sum_$filename";


$avg_1H = $avg_1M = $avg_1L = $avg_2H = $avg_2M = $avg_2L = $avg_3H = $avg_3M = $avg_3L = $avg_4H = $avg_4M = $avg_4L = $avg_5H = $avg_5M = $avg_5L = $avg_6H = $avg_6M = $avg_6L = $avg_7H = $avg_7M = $avg_7L = $avg_8H = $avg_8M = $avg_8L = $avg_9H = $avg_9M = $avg_9L = $avg_10H = $avg_10M = $avg_10L = 0;



############################################################
##################  Subroutine definitions #################
############################################################
sub double_rand{
	$max = shift || 1; 
	$iv = int rand(RAND_MAX) << RANDBITS |
	      int rand(RAND_MAX); 
	return $max * ($iv/2**(2*RANDBITS)); 
}

sub calcAlleleFreqs(\@@) { 
    # Used to calculate allele frequencies for a given Species population.
    
    # Arguments, in order:
        # Population array for which allele frequencies will be generated (e.g., @femalesC1)
        # Frequency of the cannibal allele (e.g., freq_c1)
    
    $popF = @_[0]; #Female population array 
    $popLenF = @$popF; # Female population size
     
    # Reminders:
    # @_[1] == freq_c
    
    ##The population array consists of three levels. Level 1 = nucleus (0) or cytoplasm (1). Level 2 = locus (e.g., 0 = Locus 0). Level 3 = chromosome 0 or chromosome 1. 
    ##Cytotype located in [1][0][0]. 
    ##Cannibal allele located at [0][0][0] and/or [0][0][1]. 

    
    ### Get counts for how often the cannibal allele occur
    # (using @_[#] to change the actual frequency variable provided, as opposed to extracting results later)
    foreach(@$popF) {    
        #Non-cannibal = 0/0
        #Cannibal = 1/1 or 0/1 
        
        ###Locus 0--location of cannibal allele       
       
        if($$_[0][0][0] == 1) {@_[1]++;}
        if($$_[0][0][1] == 1) {@_[1]++;}               
    }
    
    ### Convert counts to frequencies
    if($popLenF > 0) {
        @_[1] = @_[1]/(2*$popLenF);
    }  
}

############################################################
###################  Counting Subroutines ##################
############################################################


sub countHighCan(\@@){

    # This subroutine will return the number of high dispersal cannibals in a given population array 
    
    # Arguments, in order: 
        # Population array for which a cannibal count will be generated 
        
    $pop = @_[0]; 

    $high_c = 0;
    foreach(@$pop){
        if($$_[0][0][0] == 1 || $$_[0][0][1] == 1){  #print "check high_c";# the individual is a cannibal
            if($$_[0][1][0] == 1 && $$_[0][1][1] == 1) {$high_c++; } # the individual has high dispersal
        }
    }
    return $high_c;
	

}


sub countMidCan(\@@){

    # This subroutine will return the number of medium dispersal cannibals in a given population array 
    
    # Arguments, in order: 
        # Population array for which a cannibal count will be generated 
        
    $pop = @_[0]; 

    $mid_c = 0;
    foreach(@$pop){
        if($$_[0][0][0] == 1 || $$_[0][0][1] == 1){  # the individual is a cannibal
            if(($$_[0][1][0] == 1 && $$_[0][1][1] == 0) || ($$_[0][1][0] == 0 && $$_[0][1][1] == 1)) {$mid_c++;} # the individual has medium dispersal
        }
    }
    return $mid_c;

}

sub countLowCan(\@@){

    # This subroutine will return the number of low dispersal cannibals in a given population array 
    
    # Arguments, in order: 
        # Population array for which a cannibal count will be generated 
        
    $pop = @_[0]; 
    
    $low_c = 0;
    foreach(@$pop){
        if($$_[0][0][0] == 1 || $$_[0][0][1] == 1){ #print "check low_c"; # the individual is a cannibal
            if($$_[0][1][0] == 0 && $$_[0][1][1] == 0) {$low_c++;} # the individual has low dispersal
        }
    }
    return $low_c;

}


sub countHighVegan(\@@){

    # This subroutine will return the number of high dispersal cannibals in a given population array 
    
    # Arguments, in order: 
        # Population array for which a cannibal count will be generated 
        
    $pop = @_[0]; 

    $high_v = 0;
    foreach(@$pop){
        if($$_[0][0][0] == 0 && $$_[0][0][1] == 0){  #print "check high_v"; # the individual is a vegan
            if($$_[0][1][0] == 1 && $$_[0][1][1] == 1) {$high_v++;} # the individual has high dispersal
        }
    }
    return $high_v;
    
}

sub countMidVegan(\@@){

    # This subroutine will return the number of medium dispersal cannibals in a given population array 
    
    # Arguments, in order: 
        # Population array for which a cannibal count will be generated 
        
    $pop = @_[0]; 

    $mid_v = 0;
    foreach(@$pop){
        if($$_[0][0][0] == 0 && $$_[0][0][1] == 0){  # the individual is a vegan
            if(($$_[0][1][0] == 1 && $$_[0][1][1] == 0) || ($$_[0][1][0] == 0 && $$_[0][1][1] == 1)) {$mid_v++;} # the individual has medium dispersal
        }
    }
    return $mid_v;
    
}


sub countLowVegan(\@@){

    # This subroutine will return the number of high dispersal cannibals in a given population array 
    
    # Arguments, in order: 
        # Population array for which a cannibal count will be generated 
        
    $pop = @_[0]; 

    $low_v = 0;
    foreach(@$pop){
        if($$_[0][0][0] == 0 && $$_[0][0][1] == 0){  # the individual is a cannibal
            if($$_[0][1][0] == 0 && $$_[0][1][1] == 0) {$low_v++;} # the individual has low dispersal
        }
    }
    return $low_v;

}


############################################################
#################  End Counting Subroutines ################
############################################################


sub calcR{
    
    # This subroutine will return the number of offspring a particular female will produce when mating. 
    # The average number of offspring per female is 2.1. 
    # If desired, population recovery can increase the reproductive output of each female in conditions of low-density populations. 
       
       $size = $_[0];  
       $total_size = $_[1]; 
       
       
       ### Equations 1--here the resouce advantage exactly offsets the allee effect except in the case of invaders
       # resource advantage 
        $resource_advantage = $master_R + (($carrying_capacity - $total_size)/$carrying_capacity)**$Rexp;
           
       # allee cost 
       $allee_cost = (($carrying_capacity - $size)/$carrying_capacity)**$Rexp; 
       
       
       $effective_R = $resource_advantage - $allee_cost; 
       
       
       if($effective_R <= 0) {return 0; }
       
       $chance = rand(); 
       
       #note--this has an upper limit of 3 as the average R 
       if($effective_R > 2){
           $pct_chance_3off = $effective_R - 2;  
           $pct_chance_2off = 1 - $pct_chance_3off; 
           if($chance < $pct_chance_2off) { $R = 2; }
           else { $R = 3; }
       }
       elsif($effective_R > 1){
           $pct_chance_2off = $effective_R - 1; 
           $pct_chance_1off = 1 - $pct_chance_2off; 
           if($chance < $pct_chance_1off) {$R = 1; }
           else {$R = 2; }
       }
       else{
           $pct_chance_1off = $effective_R; 
           $pct_chance_0off = 1 - $pct_chance_1off;  
           if($chance < $pct_chance_0off) {$R = 0; }
           else {$R = 1; } 
       }

      
       return $R;   
}

sub calcCarryingCapacity {
    # Used to calculate maximum number of offspring produced by a given generation (maximum = if no fitness costs are applied and if cannibals don't eat all the adults before the mating season is over).
    # If applicable, the maximum number of offspring is then scaled by rare male fitness cost outside of the subroutine.
    # Final offspring number is possibly reduced again as a result of fitness costs within the createCarnivoreOffspring subroutine.
    
    ### Aka--this subroutine returns the carrying capacity of the habitat. 
    ###         The actual number of offspring produced that generation will depend on the number of females and the reproductive output of each female in that habitat. 
        
    return $carrying_capacity; 
}

sub calcSpeciesProportion {
    # Used to calculate the proportion of one species in the habitat 
    
    $one_species = $_[0]; # Number of one species in the habitat (e.g., $totalA1)
    $both_species = $_[1]; # Total number of individuals in the habitat (e.g., $total1) 
    
    if($both_species > 0) {$proportion = $one_species/$both_species;} else {$proportion = 0;} 
    
    return $proportion;
}

sub calcCannibalRatio(\@@) {
    # Used to calculate proportion of Carnivore females with cannibal genotype 
    
    # Arguments, in order:
        # Array containing genetic information of females of population    
        
    $popF = @_[0]; #Female population array 
    $popLenF = @$popF; # Female population size
    
    $cannibals = 0; 
    

    foreach(@$popF) {    
        #Non-cannibal = 0/0
        #Cannibal = 1/1 or 0/1 
        
        ###Locus 0--location of cannibal allele       
       
        if($$_[0][0][0] == 1 || $$_[0][0][1] == 1) {$cannibals++;}              
    }
    
    if($popLenF > 0) {$cannibals = $cannibals/$popLenF;} 
    return $cannibals;         
}

sub createVeganOffspring (\@\@) {
    # Creates a Vegan offspring and pushes it to the appropriate offspring array.
    
    # Arguments, in order:
        # Array that female offspring will be sent to (e.g., @femaleOffspring1).
        # Array that male offspring will be sent to (e.g., @maleOffspring1).
    
    $pushFemale = @_[0];
    $pushMale = @_[1];
    
    $gender = rand();
    if ($gender > $offSR) { push(@$pushMale, 0);} # Offspring is male
    else { push(@$pushFemale, 0); } # Offspring is female
}



sub createCarnivoreOffspring (\@\@\@\@) {
    # Creates a Carnivore offspring from two parents and pushes it to the appropriate offspring array.
    
    # Arguments, in order:
        # Array containing the mother's genetic information (e.g., @{$femalesC1[0]}).
        # Array containing the father's genetic information (e.g., @{$malesC1[0]}).
        # Array that female offspring will be sent to (e.g., @FoffspringC1).
        # Array that male offspring will be sent to (e.g., @MoffspringC1).
    
    $mother = @_[0];
    $father = @_[1];
    $pushFemale = @_[2];
    $pushMale = @_[3];
    
    # Create offspring nucleus and cytoplasm 
    @offN = ();
    @offC = (); 
  
    ##################################################
    ############ Cytoplasmic inheritance #############
    ##################################################
    
    $offC[0][0][0] = 0;
      
    ##################################################
    ############## Nuclear inheritance ###############
    ##################################################
    
    ### Locus 0, Cannibal
    # Determine which chromosome (0/1) is inherited from each parent
    $ch0F = int(rand(2)); # From mother
    $ch0M = int(rand(2)); # From father
    # Send relevant allele to offspring nucleus array
    push(@{$offN[0][0]}, $$mother[0][0][$ch0F]);
    push(@{$offN[0][0]}, $$father[0][0][$ch0M]);


    ### Locus 1, Dispersal
    # Determine which chromosome (0/1) is inherited from each parent
    $ch0F = int(rand(2)); # From mother
    $ch0M = int(rand(2)); # From father
    # Send relevant allele to offspring nucleus array
    push(@{$offN[0][1]}, $$mother[0][1][$ch0F]);
    push(@{$offN[0][1]}, $$father[0][1][$ch0M]);
    
    
    ##################################################
    #### Potential nuclear mutation for Cannibal #####
    ##################################################
       
    if($mutrt > 0) { $mutate0 = int(double_rand((1/$mutrt)));  }
    else { $mutate0 = 0; }
    if($mutrt > 0) { $mutate1 = int(double_rand((1/$mutrt)));  }
    else { $mutate1 = 0; }
        
        
    if($mutate0 == 42) { 
         $offN[0][0][0] = 1 - $offN[0][0][0];} #1->0 AND 0->1 
    if($mutate1 == 42) { 
         $offN[0][0][1] = 1 - $offN[0][0][1];} #1->0 AND 0->1 
   
        
    ##################################################
    ###### Preparation for offspring production ######
    ##################################################
    
    $gender = rand();
    $success = rand(); 
    
    ### Create indicators for whether mother has the cannibal allele
    if($$mother[0][0][0] == 1 && $$mother[0][0][1] == 1) {$maternalCan = 2;} 
    elsif($$mother[0][0][0] == 1 || $$mother[0][0][1] == 1) {$maternalCan = 1;}
    else {$maternalCan = 0;} 


    ### Create indicators for whether mother has the dispersal allele
    if($$mother[0][1][0] == 1 && $$mother[0][1][1] == 1) {$maternalDis = 2;} 
    elsif($$mother[0][1][0] == 1 || $$mother[0][1][1] == 1) {$maternalDis = 1;}
    else {$maternalDis = 0;} 

    
    ##################################################
    ############## Offspring generation ##############
    ##################################################

    # Mother is homozygous for cannibalism 
    if($maternalCan == 2 && $success > $canFC){
        if($gender > $offSR) { push(@$pushMale, [@offN, @offC]); return 1;   } #offspring is male
        else { push(@$pushFemale, [@offN, @offC]); return 2;  } #offspring is female 
    }
    
    # Mother is heterozygous for cannibalism 
    elsif($maternalCan == 1 && $success > $canFC){
        if($gender > $offSR) { push(@$pushMale, [@offN, @offC]); return 1;} #offspring is male
        else { push(@$pushFemale, [@offN, @offC]); return 2; } #offspring is female 
    }
    
    # Mother is homozygous for non-cannibalism 
    elsif($success > $flatFC){
        if($gender > $offSR) { push(@$pushMale, [@offN, @offC]); return 1; } #offspring is male
        else { push(@$pushFemale, [@offN, @offC]); return 2; } #offspring is female 
    }    
      
}

sub immigrationRates {  
    # Used in migration() to determine immigration rate from an origin habitat to a given destination habitat.  
    
    # Arguments, in order: 
        # origin (habitat number)
        # destination (habitat number)
        # dispersal code (e.g, 'H')
    
    $origin = $_[0];
    $destination = $_[1];
    $dispersal_code = $_[2]; 
    
    if($dispersal_code eq 'H') {
        if($origin == $destination) {$immRate = $staySameDisH;} 
        else {$immRate = $moveNewDisH; }
    }
    elsif($dispersal_code eq 'M') {
        if($origin == $destination) {$immRate = $staySameDisM;}
        else {$immRate = $moveNewDisM; }
    }
    else{
        if($origin == $destination) {$immRate = $staySameDisL;}
        else {$immRate = $moveNewDisL; }
    }
    
    return $immRate;
}


sub migration ($\@\@\@\@\@\@\@\@\@\@\@\@\@\@\@\@\@\@\@\@\@\@\@\@\@) {
    # Used to randomly move individuals from one habitat to another.
    
    # Arguments, in order:
        # Habitat number (e.g., 1)
 
        # Origin array (e.g., @FoffspringA1)
        # Destination array for habitat 1 (e.g., @FmigrantsA1)
        # Destination array for habitat 2 (e.g., @FmigrantsA2)
        # Destination array for habitat 3 (e.g., @FmigrantsA3)
        # Destination array for habitat 4 (e.g., @FmigrantsA4)
        # Destination array for habitat 5 (e.g., @FmigrantsA5)
        # Destination array for habitat 6 (e.g., @FmigrantsA6)
        # Destination array for habitat 7 (e.g., @FmigrantsA7)
        # Destination array for habitat 8 (e.g., @FmigrantsA8)
        # Destination array for habitat 9 (e.g., @FmigrantsA9)
        # Destination array for habitat 10 (e.g., @FmigrantsA10)
    
    $habitatNumber = $_[0];
    
    $mothers = @_[1]; 
    $female_tracker = @_[2]; 
    $male_tracker = @_[3]; 
       
    $push_FemH1 = @_[4]; 
    $push_FemH2 = @_[5]; 
    $push_FemH3 = @_[6]; 
    $push_FemH4 = @_[7]; 
    $push_FemH5 = @_[8]; 
    $push_FemH6 = @_[9]; 
    $push_FemH7 = @_[10]; 
    $push_FemH8 = @_[11]; 
    $push_FemH9 = @_[12]; 
    $push_FemH10 = @_[13]; 
    
    $push_MalH1 = @_[14]; 
    $push_MalH2 = @_[15]; 
    $push_MalH3 = @_[16]; 
    $push_MalH4 = @_[17]; 
    $push_MalH5 = @_[18]; 
    $push_MalH6 = @_[19]; 
    $push_MalH7 = @_[20]; 
    $push_MalH8 = @_[21]; 
    $push_MalH9 = @_[22]; 
    $push_MalH10 = @_[23]; 
    
    $femaleOffspring = @_[24]; 
    $maleOffspring = @_[25]; 
    
    $num_mothers = @$mothers; 
    
    
    
    $current_female = 0; $current_male = 0; 

    ### Migrate each individual to new array (which could still be native habitat) 
    for($mother = 0; $mother < $num_mothers; $mother++) { 
        $migrate = double_rand();


        # Find Mother's Dispersal Phenotype
        if($$_[0][0][0] == 1 && $$_[0][0][1] == 1) {$dispersal_code = 'H';} 
        elsif($$_[0][0][0] == 1 || $$_[0][0][1] == 1) {$dispersal_code = 'M';}
        else {$dispersal_code = 'L';}
        
        # check the tracker numbers for this mother        
        $num_female_children = $female_tracker->[$mother];  
        $num_male_children = $male_tracker->[$mother]; 
        #print "children: $num_female_children and $num_male_children\n"; 

        

        ### Based on immigration rates, set probabilites for sending an individual to each habitat
        $imm1 = immigrationRates($habitatNumber, 1, $dispersal_code);
        $imm2 = $imm1 + immigrationRates($habitatNumber, 2, $dispersal_code);
        $imm3 = $imm2 + immigrationRates($habitatNumber, 3, $dispersal_code);
        $imm4 = $imm3 + immigrationRates($habitatNumber, 4, $dispersal_code);
        $imm5 = $imm4 + immigrationRates($habitatNumber, 5, $dispersal_code);
        $imm6 = $imm5 + immigrationRates($habitatNumber, 6, $dispersal_code);
        $imm7 = $imm6 + immigrationRates($habitatNumber, 7, $dispersal_code);
        $imm8 = $imm7 + immigrationRates($habitatNumber, 8, $dispersal_code); 
        $imm9 = $imm8 + immigrationRates($habitatNumber, 9, $dispersal_code); 


        if($migrate < $imm1) {
           

            # migrate those females from the female array 
            for($move = 0; $move < $num_female_children; $move++){ 
                push(@$push_FemH1, $femaleOffspring->[$current_female]); 
                $current_female++; 
            }
            
            # migrate those males from the male array 
            for($move = 0; $move < $num_male_children; $move++){
                push(@$push_MalH1, $maleOffspring->[$current_male]); 
                $current_male++; 
            }
        }
                       
        elsif($migrate >= $imm1 && $migrate < $imm2) {
            # migrate those females from the female array 
            for($move = 0; $move < $num_female_children; $move++){ 
                push(@$push_FemH2, $femaleOffspring->[$current_female]); 
                $current_female++; 
            }
            
            # migrate those males from the male array 
            for($move = 0; $move < $num_male_children; $move++){
                push(@$push_MalH2, $maleOffspring->[$current_male]); 
                $current_male++; 
            }
        }
            
        elsif($migrate >= $imm2 && $migrate < $imm3) {
            # migrate those females from the female array 
            for($move = 0; $move < $num_female_children; $move++){ 
                push(@$push_FemH3, $femaleOffspring->[$current_female]); 
                $current_female++; 
            }
            
            # migrate those males from the male array 
            for($move = 0; $move < $num_male_children; $move++){
                push(@$push_MalH3, $maleOffspring->[$current_male]); 
                $current_male++; 
            }
        }
            
        elsif($migrate >= $imm3 && $migrate < $imm4) {
            # migrate those females from the female array 
            for($move = 0; $move < $num_female_children; $move++){ 
                push(@$push_FemH4, $femaleOffspring->[$current_female]); 
                $current_female++; 
            }
            
            # migrate those males from the male array 
            for($move = 0; $move < $num_male_children; $move++){
                push(@$push_MalH4, $maleOffspring->[$current_male]); 
                $current_male++; 
            }
        }
            
        elsif($migrate >= $imm4 && $migrate < $imm5) {
            # migrate those females from the female array 
            for($move = 0; $move < $num_female_children; $move++){ 
                push(@$push_FemH5, $femaleOffspring->[$current_female]); 
                $current_female++; 
            }
            
            # migrate those males from the male array 
            for($move = 0; $move < $num_male_children; $move++){
                push(@$push_MalH5, $maleOffspring->[$current_male]); 
                $current_male++; 
            }
        }
            
        elsif($migrate >= $imm5 && $migrate < $imm6) {
            # migrate those females from the female array 
            for($move = 0; $move < $num_female_children; $move++){ 
                push(@$push_FemH6, $femaleOffspring->[$current_female]); 
                $current_female++; 
            }
            
            # migrate those males from the male array 
            for($move = 0; $move < $num_male_children; $move++){
                push(@$push_MalH6, $maleOffspring->[$current_male]); 
                $current_male++; 
            }
        }
            
        elsif($migrate >= $imm6 && $migrate < $imm7) {
            # migrate those females from the female array 
            for($move = 0; $move < $num_female_children; $move++){ 
                push(@$push_FemH7, $femaleOffspring->[$current_female]); 
                $current_female++; 
            }
            
            # migrate those males from the male array 
            for($move = 0; $move < $num_male_children; $move++){
                push(@$push_MalH7, $maleOffspring->[$current_male]); 
                $current_male++; 
            }
        }
            
        elsif($migrate >= $imm7 && $migrate < $imm8) {
            # migrate those females from the female array 
            for($move = 0; $move < $num_female_children; $move++){ 
                push(@$push_FemH8, $femaleOffspring->[$current_female]); 
                $current_female++; 
            }
            
            # migrate those males from the male array 
            for($move = 0; $move < $num_male_children; $move++){
                push(@$push_MalH8, $maleOffspring->[$current_male]); 
                $current_male++; 
            }
        }
            
        elsif($migrate >= $imm8 && $migrate < $imm9) {
            # migrate those females from the female array 
            for($move = 0; $move < $num_female_children; $move++){ 
                push(@$push_FemH9, $femaleOffspring->[$current_female]); 
                $current_female++; 
            }
            
            # migrate those males from the male array 
            for($move = 0; $move < $num_male_children; $move++){
                push(@$push_MalH9, $maleOffspring->[$current_male]); 
                $current_male++; 
            }
        }
            
        else {
            # migrate those females from the female array 
            for($move = 0; $move < $num_female_children; $move++){ 
                push(@$push_FemH10, $femaleOffspring->[$current_female]); 
                $current_female++; 
            }
            
            # migrate those males from the male array 
            for($move = 0; $move < $num_male_children; $move++){
                push(@$push_MalH10, $maleOffspring->[$current_male]); 
                $current_male++; 
            }
        }
        
    }
}

sub shiftPopulation(\@\@\@\@) { 
    # Used to move individuals from one array into another array (e.g., to make offspring into the next generation of adults)
    
    # Arguments, in order:
        # Destination array A (e.g., @femalesC4)
        # Origin array A (e.g., @FoffspringC4)
        # Destination array B (e.g., @malesC4)
        # Origin array B (e.g., @MoffspringC4)
    
    $moveTo_A = @_[0];
    $moveFrom_A = @_[1];
    $moveTo_B = @_[2];
    $moveFrom_B = @_[3];
    
    foreach(@$moveFrom_A) {push(@$moveTo_A, $_);}
    foreach(@$moveFrom_B) {push(@$moveTo_B, $_);}
}

            
sub makeNextGenCompetingSpecies($$$$$$$\@\@\@\@\@\@\@\@\@\@\@\@) {
    
    # This subroutine goes through the procedure for a given generation's mating season. 
    # It takes care of both Species A and Species B in the same habitat with a shared carrying capacity. 
    
    # Arguments, in order: 
        # Number of High Dispersal females (e.g., $lenF1H)
        # Number of High Dispersal males (e.g., $lenM1H)
        # Number of Medium Dispersal females (e.g., $lenF1M)
        # Number of Medium Dispersal males (e.g., $lenM1M)
        # Number of Low Dispersal males (e.g., $lenF1L)
        # Number of Low Dispersal males (e.g., $lenM1L)
        # Maximum offspring (e.g., $max_offspring1)
        # Female Species A adults array (e.g., @femalesA1)
        # Male Species A adults array (e.g., @malesA1)
        # Female Species B adults array (e.g., @femalesB1)
        # Male Species B adults array (e.g., @malesB1)
        # Female Species A offspring array (e.g., @FoffspringA1)
        # Male Species A offspring array (e.g., @MoffspringA1)        
        # Female Species B offspring array (e.g., @FoffspringB1)
        # Male Species B offspring array (e.g., @MoffspringA1) 
        # female offspring tracker Species A
        # male offspring tracker Species A         
        # female offspring tracker Species B
        # male offspring tracker Species B 
           
        $lenFemH = $_[0];
        $lenMalH = $_[1];  
        $lenFemM = $_[2]; 
        $lenMalM = $_[3]; 
        $lenFemL = $_[4]; 
        $lenMalL = $_[5]; 
        $max_off = $_[6]; 
        $femalesA = @_[7]; 
        $malesA = @_[8];
        $femalesB = @_[9]; 
        $malesB = @_[10];  
        $female_offspringA = @_[11]; 
        $male_offspringA = @_[12]; 
        $female_offspringB = @_[13]; 
        $male_offspringB = @_[14]; 
        $female_trackerA = @_[15]; 
        $male_trackerA = @_[16]; 
        $female_trackerB = @_[17]; 
        $male_trackerB = @_[18]; 
        
        
        
        $females_remaining_A = $lenFemH + $lenFemM + $lenFemL;  #number of females left to produce offspring 
        $lenFemA = $females_remaining_A;
        $lenMalA = $lenMalH + $lenMalM + $lenMalL;
        $off = 0; 
        
        # mating will continue until either the carrying capcity is reached or until all females from both species have produced offspring 
        while($off < $max_off && ($females_remaining_A > 0 || $females_remaining_B > 0)){
            
            $sizeA = $lenFemA + $lenMalA; 

                if($lenFemA == 0 || $lenMalA == 0) {$females_remaining_A = 0;}
                if($females_remaining_A > 0){ #Species A can reproduce unless all their females have already produced offspring 
                    
                    # produce a Species A offspring 
                    
                    #find the female and male mates 
                    $Fmate = int(rand($lenFemA));
                    $Mmate = int(rand($lenMalA));
                    
                    #calculate reproductive rate of that female
                    $size = $lenFemA + $lenMalA;  
                    $total = $size + $lenFemB + $lenMalB; 
                    $R = calcR($size, $total); 
                    
                    $male_countA = 0; $female_countA = 0; 
                    #make offspring (based on that reproductive rate) whether or not Species B female is cannibal or non-cannibal 
                    for($off3 = 0; $off3 < $R; $off3++){
                        $track_gender = createCarnivoreOffspring(@{$$femalesA[$Fmate]}, @{$$malesA[$Mmate]}, @$female_offspringA, @$male_offspringA); 
                        $off++; 
                        if($track_gender == 1) {$male_countA++; }
                        else {$female_countA++; }
                    }
                    
                    ### Does the mother in question possess a cannibal genotype? 1=yes; 0=no
                    if($$femalesA[$Fmate][0][0][0] == 1 || $$femalesA[$Fmate][0][0][1] == 1) {$maternalCan = 1;} else {$maternalCan = 0;}
                    
                    #if female is cannibal, then she will make additional offspring 
                    $ChanceCanOff = rand(); 
                    $ChanceCanCost = rand();  
                    $cancost = 1; #$cancost = rand($cancost_rand);
                    
                    # If the mother expresses the cannibal phenotype... 
                    if($maternalCan == 1 && $ChanceCanOff < $chance_cannibalism_offspring){ #if cannibal mother produces additional offspring  
                        # Cannibal mother's additional offspring
                        for($add_off = 0; $add_off < $num_cannibal_offspring; $add_off++){ #make the additional offspring that a cannibal produces
                            $track_gender = createCarnivoreOffspring(@{$$femalesA[$Fmate]}, @{$$malesA[$Mmate]}, @$female_offspringA, @$male_offspringA); #Second Offspring if the above is true 
                            $off++;
                            if($track_gender == 1) {$male_countA++; }
                            else {$female_countA++; }  
                        }
                    }     
                    
                    if($maternalCan == 1 && $ChanceCanCost < $chance_cannibalism_cost){ #if cannibal mother kills off others of her population 
                        # Individuals killed by cannibal mother when gathering resources to lay eggs 
                        for($die = 0; $die < $cancost; $die++){  
                            $gender = int(rand(2));
                            if($gender == 0){$lenFemA = $lenFemA - 1; $females_remaining_A -= 1; #she can kill off either males or females, and they are no longer able to mate (because they are dead) 
                                }else{$lenMalA = $lenMalA - 1;}
                        } 
                    }
                    if($lenMalA <= 0 || $lenFemA <= 0){$females_remaining_A = 0; }                    
                    
                    $females_remaining_A -= 1; 
                    push(@$female_trackerA, $female_countA); 
                    push(@$male_trackerA, $male_countA); 
                }

        }
}    
 



############################################################
####################### File Setup #########################
############################################################

### Create files to save output

    # Simulation summary (multiple trials)- used with parameter loops
    open(FHd, ">$sumname.csv");
	open(FH, ">G:/Shared drives/sim1119/onespecies/$sumname.csv");

	select((select(FH), $|=1)[0]);
	open(FH5d, ">$sumname.columns.csv");
	open(FH5, ">G:/Shared drives/sim1119/onespecies/$sumname.columns.csv");

	select((select(FH5), $|=1)[0]);
    # Per-simulation results
    #open(FH1, ">scn5v3_eachSim.csv");
    # Per generation results
    #open(FH2, ">scn5v3_gens.csv"); 


### Print non-looping parameters and column descriptions to files
# Simulation summary results (use skip __ lines when reading in file, e.g. to R)
print(FH "### PARAMETERS ###\nScenario number = 5 (interspecific competition/metapopulation/migration)\nBeginning of migration = generation $invgen\nEnding generation = $endgen\nReproductive output = $femrep (control by: 1=females / 0=both sexes)\nOffspring sex ratio = $offSR\nFC scaling = $flatFC/$canFC0/$canFC1/$canFC2 (scalars)\n\ncarrying capacity = $carrying_capacity\nmutation rate = $mutrt\nmigration rate H = $migrtH\nmigration rate M = $migrtM\nmigration rate L = $migrtL\ncancost = $cancost_rand\ncan_off = $num_cannibal_offspring\ndead = $dead\n\n");
print(FHd "### PARAMETERS ###\nScenario number = 5 (interspecific competition/metapopulation/migration)\nBeginning of migration = generation $invgen\nEnding generation = $endgen\nReproductive output = $femrep (control by: 1=females / 0=both sexes)\nOffspring sex ratio = $offSR\nFC scaling = $flatFC/$canFC0/$canFC1/$canFC2 (scalars)\n\ncarrying capacity = $carrying_capacity\nmutation rate = $mutrt\nmigration rate H = $migrtH\nmigration rate M = $migrtM\nmigration rate L = $migrtL\ncancost = $cancost_rand\ncan_off = $num_cannibal_offspring\ndead = $dead\n\n");

#print(FH "### COLUMNS ###\navg_total = Average final individuals\navg_totalA = Average final Species A individuals\navg_totalB = Average final Species B individuals\navgA/B# = Average number of Species A or B in Habitat #\navg_propA = average proportion of Species A in metacommunity\navg_propB = average proportion of Species B in metacommunity\navg_endgen = the average ending generation among runs\nruns = Number of simulations run with given parameter set\n\n### END HEADER ###\n\n");

# Per-simulation results (skip=_)
#print(FH1 "### PARAMETERS ###\nScenario number = 5 (interspecific competition/metapopulation/migration)\nBeginning of migration = generation $invgen\nEnding generation = $endgen\nReproductive output = $femrep (control by: 1=females / 0=both sexes)\nOffspring sex ratio = $offSR\nFC scaling = $flatFC/$canFC0/$canFC1/$canFC2 (scalars)\n\ncarrying capacity = $carrying_capacity\nmutation rate = $mutrt\nmigration rate A = $migrtA\nmigration rate B = $migrtB\ncancost = $cancost_rand\ncan_off = $num_cannibal_offspring\ndead = $dead\n\n");
#print(FH1 "### COLUMNS ###\ntotal = # individuals when simulation ends\ntotalA = # of Species A when simulation ends\ntotalB = # of Species B when simulation ends\ntotalA# = number of Species A in Habitat #\ntotalB# = number of Species B in Habitat #\npropA = proportion of Species A in metacommunity\npropB = proportion of Species B in metacommunity\nendgen = Final generation\nrun = Run number with that set of parameters\n\n### END HEADER ###\n\n");

# Per generation results (skip=_)
#print(FH2 "### PARAMETERS ###\nScenario number = 5 (interspecific competition/metapopulation/migration)\nBeginning of migration = generation $invgen\nEnding generation = $endgen\nReproductive output = $femrep (control by: 1=females / 0=both sexes)\nOffspring sex ratio = $offSR\nFC scaling = $flatFC/$canFC0/$canFC1/$canFC2 (scalars)\n\ncarrying capacity = $carrying_capacity\nmutation rate = $mutrt\nmigration rate A = $migrtA\nmigration rate B = $migrtB\ncancost = $cancost_rand\ncan_off = $num_cannibal_offspring\ndead = $dead\n\n");
#print(FH2 "### COLUMNS ###\nTotalA/B#c/n = number of cannibal/non-cannibal individuals of Species A or B in Habitat\navgA/B# = avg number of individuals of Species A or B in Habitat #\ngen = generation #\nrun = Run number of the simulation\n\n### END HEADER ###\n\n");

# Per generation results (skip=_)
#print(FH3 "### PARAMETERS ###\nScenario number = 5 (interspecific competition/metapopulation/migration)\nBeginning of migration = generation $invgen\nEnding generation = $endgen\nReproductive output = $femrep (control by: 1=females / 0=both sexes)\nOffspring sex ratio = $offSR\nFC scaling = $flatFC/$canFC0/$canFC1/$canFC2 (scalars)\n\ncarrying capacity = $carrying_capacity\nmutation rate = $mutrt\nmigration rate A = $migrtA\nmigration rate B = $migrtB\ncancost = $cancost_rand\ncan_off = $num_cannibal_offspring\ndead = $dead\n\n");
#print(FH3 "### COLUMNS ###\nTotalA/B#c/n = number of cannibal/non-cannibal individuals of Species A or B in Habitat\navgA/B# = avg number of individuals of Species A or B in Habitat #\ngen = generation #\nrun = Run number of the simulation\n\n### END HEADER ###\n\n");


### Print column name headers to files
#print(FH "simno,end_A1H,end_A1M,end_A1L,end_A2H,end_A2M,end_A2L,end_A3H,end_A3M,end_A3L,end_A4H,end_A4M,end_A4L,end_A5H,end_A5M,end_A5L,end_A6H,end_A6M,end_A6L,end_A7H,end_A7M,end_A7L,end_A8H,end_A8M,end_A8L,end_A9H,end_A9M,end_A9L,end_A10H,end_A10M,end_A10L\n");



#print(FH1 "total,totalA,totalB,totalA1,totalA2,totalA3,totalA4,totalA5,totalA6,totalA7,totalA8,totalA9,totalA10,totalB1,totalB2,totalB3,totalB4,totalB5,totalB6,totalB7,totalB8,totalB9,totalB10,propA,propB,endgen,run\n");
#print(FH2 "totalA1c,totalA1n,totalA2c,totalA2n,totalA3c,totalA3n,totalA4c,totalA4n,totalA5c,totalA5n,totalA6c,totalA6n,totalA7c,totalA7n,totalA8c,totalA8n,totalA9c,totalA9n,totalA10c,totalA10n,totalB1c,totalB1n,totalB2c,totalB2n,totalB3c,totalB3n,totalB4c,totalB4n,totalB5c,totalB5n,totalB6c,totalB6n,totalB7c,totalB7n,totalB8c,totalB8n,totalB9c,totalB9n,totalB10c,totalB10n,gen,run\n"); 
#print(FH3 "totalA1c,totalA1n,totalA2c,totalA2n,totalA3c,totalA3n,totalA4c,totalA4n,totalA5c,totalA5n,totalA6c,totalA6n,totalA7c,totalA7n,totalA8c,totalA8n,totalA9c,totalA9n,totalA10c,totalA10n,totalB1c,totalB1n,totalB2c,totalB2n,totalB3c,totalB3n,totalB4c,totalB4n,totalB5c,totalB5n,totalB6c,totalB6n,totalB7c,totalB7n,totalB8c,totalB8n,totalB9c,totalB9n,totalB10c,totalB10n,gen,run\n"); 





####################### Simulations ########################
############################################################

##### OVERVIEW: Every simulation... #####
    # 1. Create initial individuals.
    # 2. Determine allele/cytoplasm frequencies.
    # 3. Assign randomly-selected females to a random mate and (possibly) create offspring from that pair.
    # 4. Put each new offspring into a sex-specific array.
    # 5. Allow individuals to migrate to new habitats. 

### Set up parameter ranges for loops, if applicable 



### Insert parameter loops, if applicable  
### File looping 



### This loop includes parameter loops

for($loop_cost = 0; $loop_cost<1; $loop_cost++){

    #if($loop_cost == 0)  {open(FH2, ">capacity90.csv"); }
    #if($loop_cost == 1)  {open(FH2, ">capacity100.csv"); }
    
    #print(FH2 "totalA1c,totalA1n,totalA2c,totalA2n,totalA3c,totalA3n,totalA4c,totalA4n,totalA5c,totalA5n,totalA6c,totalA6n,totalA7c,totalA7n,totalA8c,totalA8n,totalA9c,totalA9n,totalA10c,totalA10n,totalB1c,totalB1n,totalB2c,totalB2n,totalB3c,totalB3n,totalB4c,totalB4n,totalB5c,totalB5n,totalB6c,totalB6n,totalB7c,totalB7n,totalB8c,totalB8n,totalB9c,totalB9n,totalB10c,totalB10n,cost,benefit,gen,run\n"); 


for($loop_benefit = 0; $loop_benefit<1; $loop_benefit++){     
    
        
    $chance_cannibalism_cost = $all_loop_costs[$loop_cost];  #probability that a cannibal mother will eat others of her population 
    $chance_cannibalism_offspring = $all_loop_benefits[$loop_benefit];  #probability that a cannibal mother will have additional offspring 


### New Data Sheet by Caroline

open(FH4d, ">$filename.csv");
open(FH4, ">G:/Shared drives/sim1119/onespecies/$filename.csv");

select((select(FH4), $|=1)[0]);
print(FH4 "gen,run,cost,benefit,migrtA,migrtB,habitat,dispersal,cannibal,pop\n");    
print(FH4d "gen,run,cost,benefit,migrtA,migrtB,habitat,dispersal,cannibal,pop\n");    
    
    
### Population survival indicators
$extinctionA = $extinctionB = 0; # simulation ends due to extinction of either Species A or Species B 
$total_extinctionsA = $total_extinctionsB = 0; # number of times Species A or Species B went extinct 
$ct=0;
$simno = 1;

$avg_total = $avg_total_A = $avg_total_B = $avg_prop_A = $avg_prop_B = $avg_endgen = 0; 
$avgA1 = $avgA2 = $avgA3 = $avgA4 = $avgA5 = $avgA6 = $avgA7 = $avgA8 = $avgA9 = $avgB1 = $avgB2 = $avgB3 = $avgB4 = $avgB5 = $avgB6 = $avgB7 = $avgB8 = $avgB9 = 0; 
$master_freq_H = $master_freq_M = $master_freq_L = 0;

# Real Terminal Output
### Run through desired number of simulations
while ($simno <= $runs) {
    #print "File: $file, Simulation $simno:\n"; 
    print "Simulation $simno:\tMutation $mutrt; Migration H $migrtH; Migration M $migrtM; Migration L $migrtL; Cost $chance_cannibalism_cost; Offspring $chance_cannibalism_offspring; R exponent $Rexp; Capacity $carrying_capacity; File $filename\n"; 
  

#print(FH5"$simno,$propA1c,$end_1\n,$propA2c,$end_2\n,$propA3c,$end_3\n,$propA4c,$end_4\n,$propA5c,$end_5\n,$propA6c,$end_6\n,$propA7c,$end_7\n,$propA8c,$end_8\n,$propA9c,$end_9\n,$propA10c,$end_10\n");


# # # Male Cannibals? Lets find out # # #    
    #if($run/25 == %0) {print "Male Cannibals at Run $runs Species A\n
   	#Habitat 1: $lenMA1c\tHabitat 2: $lenMA2c\tHabitat 3: $lenMA3c\tHabitat 4: $lenMA4c\tHabitat 5: $lenMA5c\tHabitat 6: $lenMA6c\tHabitat 7: $lenMA7c\tHabitat 8: $lenMA8c\tHabitat 9: $lenMA9c\tHabitat 10: $lenMA10c\n
   	#Species B\n
   	#Habitat 1: $lenMB1c\tHabitat 2: $lenMB2c\tHabitat 3: $lenMB3c\tHabitat 4: $lenMB4c\tHabitat 5: $lenMB5c\tHabitat 6: $lenMB6c\tHabitat 7: $lenMB7c\tHabitat 8: $lenMB8c\tHabitat 9: $lenMB9c\tHabitat 10: $lenMB10c\n"}
    




    ##################################################
    ############ Set up initial populations ##########
    ##################################################
    
    $gen = 0;


    ### Create arrays for each sex
    # Name key:
        # Format: (fe/males)(A/B)(#)
        # fe/male = sex of individuals
        # A/B = Species A/Species B individuals
        # # = habitat number
    

    
    # Male and female Species A populations of Habitats 1-9 
    @femalesA1 = (); @malesA1 = (); 
    @femalesA2 = (); @malesA2 = (); 
    @femalesA3 = (); @malesA3 = (); 
    @femalesA4 = (); @malesA4 = (); 
    @femalesA5 = (); @malesA5 = (); 
    @femalesA6 = (); @malesA6 = (); 
    @femalesA7 = (); @malesA7 = (); 
    @femalesA8 = (); @malesA8 = (); 
    @femalesA9 = (); @malesA9 = (); 
    @femalesA10 = (); @malesA10 = (); 
    
    
    ##### Produce initial individuals #####
    ### Determine how many females to create
    $k = $carrying_capacity; 
    $out = $offSR * $k;
    # Ensure the output does not end exactly in 0.5, which causes odd population sizes due to rounding error
    if(($out - 0.5) =~ /\./) {$out = $out;}
    else {$out = $out + (.00001 - rand(.00002));}
    
    
    ### Reminder--how genetic information is stored 
   
    # Two-element list for each individual, contained within the population array 
        # Format:
            # [ [nucleus], [cytoplasm] ]
            # Within nucleus:       [locus 0], [locus 1]
            # > Within each locus:  allele on chromosome A, allele on chromosome B
            # Within cytoplasm:     [[cytotype]]
    
    ### Each Habitat has a carrying capacity of 100
    ### High Dispersal Phenotype: Habitas 1-5
    ### Low Dispersal Phenotype: Habitats 6-10
    
    
    # Habitat 1
    for($ind = 0; $ind < $out; $ind++){ #out=0.5*k
        if($mutrt > 0){
            $mutate0 = int(double_rand((1/$mutrt))); 
            $mutate1 = int(double_rand((1/$mutrt)));
        } 
        else { $mutate0 = 0; $mutate1 = 0; }
        if($mutate0 == 42 && $mutate1 == 42){#both canibal copies get a mutation at the same time
            push(@femalesA1, [ [[1, 1], [1, 1]], [[0]] ]); }
        elsif($mutate0 == 42 || $mutate1 == 42){
            push(@femalesA1, [ [[1, 0], [1, 1]], [[0]] ]); print "canibal mutation";}  #locus 0-i.e. canibal gets one mutation          
        else{
            push(@femalesA1, [ [[0, 0], [1, 1]], [[0]] ]); }  #no mutation for locus 0
    }
    for($ind = 0; $ind < $k-$out; $ind++){ #$k-$out=0.5*k..this also mutates the 0 locus, but for males
        if($mutrt > 0){
            $mutate0 = int(double_rand((1/$mutrt))); 
            $mutate1 = int(double_rand((1/$mutrt)));
        } 
        else { $mutate0 = 0; $mutate1 = 0; }
        if($mutate0 == 42 && $mutate1 == 42){
            push(@malesA1, [ [[1, 1], [1, 1]], [[0]] ]); }
        elsif($mutate0 == 42 || $mutate1 == 42){
            push(@malesA1, [ [[1, 0], [1, 1]], [[0]] ]); }        
        else{
            push(@malesA1, [ [[0, 0], [1, 1]], [[0]] ]); } 
    }

    ### Habitat 2
    for($ind = 0; $ind < $out; $ind++){ 
        if($mutrt > 0){
            $mutate0 = int(double_rand((1/$mutrt))); 
            $mutate1 = int(double_rand((1/$mutrt)));
        } 
        else { $mutate0 = 0; $mutate1 = 0; } 
        if($mutate0 == 42 && $mutate1 == 42){
            push(@femalesA2, [ [[1, 1], [0, 0]], [[0]] ]); }
        elsif($mutate0 == 42 || $mutate1 == 42){
            push(@femalesA2, [ [[1, 0], [0, 0]], [[0]] ]); }            
        else{
            push(@femalesA2, [ [[0, 0], [0, 0]], [[0]] ]); } 
    }
    for($ind = 0; $ind < $k-$out; $ind++){ 
        if($mutrt > 0){
            $mutate0 = int(double_rand((1/$mutrt))); 
            $mutate1 = int(double_rand((1/$mutrt)));
        } 
        else { $mutate0 = 0; $mutate1 = 0; }
        if($mutate0 == 42 && $mutate1 == 42){
            push(@malesA2, [ [[1, 1], [0, 0]], [[0]] ]); }
        elsif($mutate0 == 42 || $mutate1 == 42){
            push(@malesA2, [ [[1, 0], [0, 0]], [[0]] ]); }        
        else{
            push(@malesA2, [ [[0, 0], [0, 0]], [[0]] ]); } 
    }

    ### Habitat 3
    for($ind = 0; $ind < $out; $ind++){

        if($mutrt > 0){
            $mutate0 = int(double_rand((1/$mutrt))); 
            $mutate1 = int(double_rand((1/$mutrt)));
        } 
        else { $mutate0 = 0; $mutate1 = 0; }
        if($mutate0 == 42 && $mutate1 == 42){
            push(@femalesA3, [ [[1, 1], [1, 1]], [[0]] ]); }
        elsif($mutate0 == 42 || $mutate1 == 42){
            push(@femalesA3, [ [[1, 0], [1, 1]], [[0]] ]); }            
        else{
            push(@femalesA3, [ [[0, 0], [1, 1]], [[0]] ]); } 
    }
    for($ind = 0; $ind < $k-$out; $ind++){ 
        if($mutrt > 0){
            $mutate0 = int(double_rand((1/$mutrt))); 
            $mutate1 = int(double_rand((1/$mutrt)));
        } 
        else { $mutate0 = 0; $mutate1 = 0; }
        if($mutate0 == 42 && $mutate1 == 42){
            push(@malesA3, [ [[1, 1], [1, 1]], [[0]] ]); }
        elsif($mutate0 == 42 || $mutate1 == 42){
            push(@malesA3, [ [[1, 0], [1, 1]], [[0]] ]); }        
        else{
            push(@malesA3, [ [[0, 0], [1, 1]], [[0]] ]); } 
    }

    ### Habitat 4
    for($ind = 0; $ind < $out; $ind++){ 

        if($mutrt > 0){
            $mutate0 = int(double_rand((1/$mutrt))); 
            $mutate1 = int(double_rand((1/$mutrt)));
        } 
        else { $mutate0 = 0; $mutate1 = 0; }
        if($mutate0 == 42 && $mutate1 == 42){
            push(@femalesA4, [ [[1, 1], [1, 1]], [[0]] ]); }
        elsif($mutate0 == 42 || $mutate1 == 42){
            push(@femalesA4, [ [[1, 0], [1, 1]], [[0]] ]); }            
        else{
            push(@femalesA4, [ [[0, 0], [1, 1]], [[0]] ]); } 
    }
    for($ind = 0; $ind < $k-$out; $ind++){ 
        if($mutrt > 0){
            $mutate0 = int(double_rand((1/$mutrt))); 
            $mutate1 = int(double_rand((1/$mutrt)));
        } 
        else { $mutate0 = 0; $mutate1 = 0; }
        if($mutate0 == 42 && $mutate1 == 42){
            push(@malesA4, [ [[1, 1], [1, 1]], [[0]] ]); }
        elsif($mutate0 == 42 || $mutate1 == 42){
            push(@malesA4, [ [[1, 0], [1, 1]], [[0]] ]); }        
        else{
            push(@malesA4, [ [[0, 0], [1, 1]], [[0]] ]); } 
    }

    ### Habitat 5
    for($ind = 0; $ind < $out; $ind++){ 

        if($mutrt > 0){
            $mutate0 = int(double_rand((1/$mutrt))); 
            $mutate1 = int(double_rand((1/$mutrt)));
        } 
        else { $mutate0 = 0; $mutate1 = 0; } 
        if($mutate0 == 42 && $mutate1 == 42){
            push(@femalesA5, [ [[1, 1], [1, 1]], [[0]] ]); }
        elsif($mutate0 == 42 || $mutate1 == 42){
            push(@femalesA5, [ [[1, 0], [1, 1]], [[0]] ]); }            
        else{
            push(@femalesA5, [ [[0, 0], [1, 1]], [[0]] ]); } 
    }
    for($ind = 0; $ind < $k-$out; $ind++){ 
        if($mutrt > 0){
            $mutate0 = int(double_rand((1/$mutrt))); 
            $mutate1 = int(double_rand((1/$mutrt)));
        } 
        else { $mutate0 = 0; $mutate1 = 0; }
        if($mutate0 == 42 && $mutate1 == 42){
            push(@malesA5, [ [[1, 1], [1, 1]], [[0]] ]); }
        elsif($mutate0 == 42 || $mutate1 == 42){
            push(@malesA5, [ [[1, 0], [1, 1]], [[0]] ]); }        
        else{
            push(@malesA5, [ [[0, 0], [1, 1]], [[0]] ]); } 
    }

    # ## Habitat 6
    for($ind = 0; $ind < $out; $ind++){ 

        if($mutrt > 0){
            $mutate0 = int(double_rand((1/$mutrt))); 
            $mutate1 = int(double_rand((1/$mutrt)));
        } 
        else { $mutate0 = 0; $mutate1 = 0; } 
        if($mutate0 == 42 && $mutate1 == 42){
            push(@femalesA6, [ [[1, 1], [0, 0]], [[0]] ]); }
        elsif($mutate0 == 42 || $mutate1 == 42){
            push(@femalesA6, [ [[1, 0], [0, 0]], [[0]] ]); }            
        else{
            push(@femalesA6, [ [[0, 0], [0, 0]], [[0]] ]); } 
    }
    for($ind = 0; $ind < $k-$out; $ind++){ 
        if($mutrt > 0){
            $mutate0 = int(double_rand((1/$mutrt))); 
            $mutate1 = int(double_rand((1/$mutrt)));
        } 
        else { $mutate0 = 0; $mutate1 = 0; }
        if($mutate0 == 42 && $mutate1 == 42){
            push(@malesA6, [ [[1, 1], [0, 0]], [[0]] ]); }
        elsif($mutate0 == 42 || $mutate1 == 42){
            push(@malesA6, [ [[1, 0], [0, 0]], [[0]] ]); }        
        else{
            push(@malesA6, [ [[0, 0], [0, 0]], [[0]] ]); } 
    }

    # Habitat 7
    for($ind = 0; $ind < $out; $ind++){
        if($mutrt > 0){
            $mutate0 = int(double_rand((1/$mutrt))); 
            $mutate1 = int(double_rand((1/$mutrt)));
        } 
        else { $mutate0 = 0; $mutate1 = 0; }
        if($mutate0 == 42 && $mutate1 == 42){
            push(@femalesA7, [ [[1, 1], [1, 1]], [[0]] ]); }
        elsif($mutate0 == 42 || $mutate1 == 42){
            push(@femalesA7, [ [[1, 0], [1, 1]], [[0]] ]); }            
        else{
            push(@femalesA7, [ [[0, 0], [1, 1]], [[0]] ]); } 
    }
    for($ind = 0; $ind < $k-$out; $ind++){ 
        if($mutrt > 0){
            $mutate0 = int(double_rand((1/$mutrt))); 
            $mutate1 = int(double_rand((1/$mutrt)));
        } 
        else { $mutate0 = 0; $mutate1 = 0; }
        if($mutate0 == 42 && $mutate1 == 42){
            push(@malesA7, [ [[1, 1], [1, 1]], [[0]] ]); }
        elsif($mutate0 == 42 || $mutate1 == 42){
            push(@malesA7, [ [[1, 0], [1, 1]], [[0]] ]); }        
        else{
            push(@malesA7, [ [[0, 0], [1, 1]], [[0]] ]); } 
    }

    ### Habitat 8
    for($ind = 0; $ind < $out; $ind++){ 
		#print "check ind $ind";

        if($mutrt > 0){
            $mutate0 = int(double_rand((1/$mutrt))); 
            $mutate1 = int(double_rand((1/$mutrt)));
        } 
        else { $mutate0 = 0; $mutate1 = 0; }
        if($mutate0 == 42 && $mutate1 == 42){
            push(@femalesA8, [ [[1, 1], [0, 0]], [[0]] ]); }
        elsif($mutate0 == 42 || $mutate1 == 42){
            push(@femalesA8, [ [[1, 0], [0, 0]], [[0]] ]); }            
        else{
            push(@femalesA8, [ [[0, 0], [0, 0]], [[0]] ]); } 
    }
    for($ind = 0; $ind < $k-$out; $ind++){ 

        if($mutrt > 0){
            $mutate0 = int(double_rand((1/$mutrt))); 
            $mutate1 = int(double_rand((1/$mutrt)));
        } 
        else { $mutate0 = 0; $mutate1 = 0; }
        if($mutate0 == 42 && $mutate1 == 42){
            push(@malesA8, [ [[1, 1], [0, 0]], [[0]] ]); }
        elsif($mutate0 == 42 || $mutate1 == 42){
            push(@malesA8, [ [[1, 0], [0, 0]], [[0]] ]); }        
        else{
            push(@malesA8, [ [[0, 0], [0, 0]], [[0]] ]); } 
    }

    ### Habitat 9
    for($ind = 0; $ind < $out; $ind++){

        if($mutrt > 0){
            $mutate0 = int(double_rand((1/$mutrt))); 
            $mutate1 = int(double_rand((1/$mutrt)));
        } 
        else { $mutate0 = 0; $mutate1 = 0; }

        if($mutate0 == 42 && $mutate1 == 42){
            push(@femalesA9, [ [[1, 1], [0, 0]], [[0]] ]); }
        elsif($mutate0 == 42 || $mutate1 == 42){
            push(@femalesA9, [ [[1, 0], [0, 0]], [[0]] ]); }            
        else{
            push(@femalesA9, [ [[0, 0], [0, 0]], [[0]] ]); } 
    }
    for($ind = 0; $ind < $k-$out; $ind++){ 
        if($mutrt > 0){
            $mutate0 = int(double_rand((1/$mutrt))); 
            $mutate1 = int(double_rand((1/$mutrt)));
        } 
        else { $mutate0 = 0; $mutate1 = 0; }
        if($mutate0 == 42 && $mutate1 == 42){
            push(@malesA9, [ [[1, 1], [0, 0]], [[0]] ]); }
        elsif($mutate0 == 42 || $mutate1 == 42){
            push(@malesA9, [ [[1, 0], [0, 0]], [[0]] ]); }        
        else{
            push(@malesA9, [ [[0, 0], [0, 0]], [[0]] ]); } 
    }

    ### Habitat 10
    for($ind = 0; $ind < $out; $ind++){

        if($mutrt > 0){
            $mutate0 = int(double_rand((1/$mutrt))); 
            $mutate1 = int(double_rand((1/$mutrt)));
        } 
        else { $mutate0 = 0; $mutate1 = 0; }
        if($mutate0 == 42 && $mutate1 == 42){
            push(@femalesA10, [ [[1, 1], [0, 0]], [[0]] ]); }
        elsif($mutate0 == 42 || $mutate1 == 42){
            push(@femalesA10, [ [[1, 0], [0, 0]], [[0]] ]); }            
        else{
            push(@femalesA10, [ [[0, 0], [0, 0]], [[0]] ]); } 
    }
    for($ind = 0; $ind < $k-$out; $ind++){ 
        if($mutrt > 0){
            $mutate0 = int(double_rand((1/$mutrt))); 
			#print "heres mutate0 $mutate0";
            $mutate1 = int(double_rand((1/$mutrt)));
        } 
        else { $mutate0 = 0; $mutate1 = 0; }
        if($mutate0 == 42 && $mutate1 == 42){
            push(@malesA10, [ [[1, 1], [0, 0]], [[0]] ]); }
        elsif($mutate0 == 42 || $mutate1 == 42){
            push(@malesA10, [ [[1, 0], [0, 0]], [[0]] ]); }        
        else{
            push(@malesA10, [ [[0, 0], [0, 0]], [[0]] ]); } 
    }

     
    
    ##################################################
    ############# Run through generations ############
    ##################################################

    
    # will continue the simulation until one of the species goes extinct or until endgen 
    do {

#print "lenf1Hc=$lenF1Hc\n\n";
#print "\tGen: $gen\n";
#print "propa1c: $propA1c\n";
#print "total1: $total1\n";

    if($gen == 2000) {print "\tGen 2000\n"; }
    if($gen == 4000) {print "\tGen 4000\n"; }
    if($gen == 6000) {print "\tGen 6000\n"; }
    if($gen == 8000) {print "\tGen 8000\n"; }

        
        ### Create scalars for number of individuals of each sex in Species A and B  
        # Name key:
            # Format: len(F/M)(A/B)(#)
            # F/M = fe/male
            # A/B = Carnivore/Vegan
            # # = habitat number

        
        # Number of males and females of Species A in each Habitat 
        $lenFA1 = @femalesA1; $lenMA1 = @malesA1;
        $lenFA2 = @femalesA2; $lenMA2 = @malesA2; 
        $lenFA3 = @femalesA3; $lenMA3 = @malesA3; 
        $lenFA4 = @femalesA4; $lenMA4 = @malesA4; 
        $lenFA5 = @femalesA5; $lenMA5 = @malesA5; 
        $lenFA6 = @femalesA6; $lenMA6 = @malesA6; 
        $lenFA7 = @femalesA7; $lenMA7 = @malesA7;
        $lenFA8 = @femalesA8; $lenMA8 = @malesA8; 
        $lenFA9 = @femalesA9; $lenMA9 = @malesA9; 
        $lenFA10 = @femalesA10; $lenMA10 = @malesA10; 
        
        # Number of males and females of Species B in each Habitat
        $lenFB1 = @femalesB1; $lenMB1 = @malesB1; 
        $lenFB2 = @femalesB2; $lenMB2 = @malesB2; 
        $lenFB3 = @femalesB3; $lenMB3 = @malesB3; 
        $lenFB4 = @femalesB4; $lenMB4 = @malesB4; 
        $lenFB5 = @femalesB5; $lenMB5 = @malesB5; 
        $lenFB6 = @femalesB6; $lenMB6 = @malesB6; 
        $lenFB7 = @femalesB7; $lenMB7 = @malesB7; 
        $lenFB8 = @femalesB8; $lenMB8 = @malesB8; 
        $lenFB9 = @femalesB9; $lenMB9 = @malesB9; 
        $lenFB10 = @femalesB10; $lenMB10 = @malesB10; 


        # Frequencies of the cannibal allele in the Species A and B populations 
 
        
                
        ### Total population sizes 
        $totalA1 = $lenFA1 + $lenMA1; 
        $totalB1 = $lenFB1 + $lenMB1; 
        $total1 = $totalA1 + $totalB1;       
        $totalF1 = $lenFA1 + $lenFB1;   
        $totalM1 = $lenMA1 + $lenMB1;   
         
        $totalA2 = $lenFA2 + $lenMA2; 
        $totalB2 = $lenFB2 + $lenMB2; 
        $total2 = $totalA2 + $totalB2;        
        $totalF2 = $lenFA2 + $lenFB2;
        $totalM2 = $lenMA2 + $lenMB2;

        $totalA3 = $lenFA3 + $lenMA3; 
        $totalB3 = $lenFB3 + $lenMB3; 
        $total3 = $totalA3 + $totalB3; 
        $totalF3 = $lenFA3 + $lenFB3;
        $totalM3 = $lenMA3 + $lenMB3;

        $totalA4 = $lenFA4 + $lenMA4; 
        $totalB4 = $lenFB4 + $lenMB4;
        $total4 = $totalA4 + $totalB4;  
        $totalF4 = $lenFA4 + $lenFB4;
        $totalM4 = $lenMA4 + $lenMB4;

        $totalA5 = $lenFA5 + $lenMA5; 
        $totalB5 = $lenFB5 + $lenMB5;
        $total5 = $totalA5 + $totalB5;  
        $totalF5 = $lenFA5 + $lenFB5; 
        $totalM5 = $lenMA5 + $lenMB5; 

        $totalA6 = $lenFA6 + $lenMA6; 
        $totalB6 = $lenFB6 + $lenMB6; 
        $total6 = $totalA6 + $totalB6; 
        $totalF6 = $lenFA6 + $lenFB6; 
        $totalM6 = $lenMA6 + $lenMB6; 

        $totalA7 = $lenFA7 + $lenMA7; 
        $totalB7 = $lenFB7 + $lenMB7; 
        $total7 = $totalA7 + $totalB7; 
        $totalF7 = $lenFA7 + $lenFB7;
        $totalM7 = $lenMA7 + $lenMB7;

        $totalA8 = $lenFA8 + $lenMA8; 
        $totalB8 = $lenFB8 + $lenMB8; 
        $total8 = $totalA8 + $totalB8; 
        $totalF8 = $lenFA8 + $lenFB8; 
        $totalM8 = $lenMA8 + $lenMB8; 

        $totalA9 = $lenFA9 + $lenMA9; 
        $totalB9 = $lenFB9 + $lenMB9; 
        $total9 = $totalA9 + $totalB9; 
        $totalF9 = $lenFA9 + $lenFB9; 
        $totalM9 = $lenMA9 + $lenMB9; 

        $totalA10 = $lenFA10 + $lenMA10; 
        $totalB10 = $lenFB10 + $lenMB10; 
        $total10 = $totalA10 + $totalB10; 
        $totalF10 = $lenFA10 + $lenFB10; 
        $totalM10 = $lenMA10 + $lenMB10; 


        #Over species counts to be used in overall Carnivore Ratio calculation
        $totalAFemale = $lenFA1 + $lenFA2 + $lenFA3 + $lenFA4 + $lenFA5 + $lenFA6 + $lenFA7 + $lenFA8 + $lenFA9 + $lenFA10;
        $totalBFemale = $lenFB1 + $lenFB2 + $lenFB3 + $lenFB4 + $lenFB5 + $lenFB6 + $lenFB7 + $lenFB8 + $lenFB9 + $lenFB10;
        $totalFemale = $totalAFemale + $totalBFemale; 
        
        $totalAMale = $lenMA1 + $lenMA2 + $lenMA3 + $lenMA4 + $lenMA5 + $lenMA6 + $lenMA7 + $lenMA8 + $lenMA9 + $lenMA10;
        $totalBMale = $lenMB1 + $lenMB2 + $lenMB3 + $lenMB4 + $lenMB5 + $lenMB6 + $lenMB7 + $lenMB8 + $lenMB9 + $lenMB10;
        $totalMale = $totalAMale + $totalBMale; 
                

        #Overall speices counts        
        $totalA = $totalAFemale + $totalAMale; 
        $totalB = $totalBFemale + $totalBMale; 
        $total = $total1 + $total2 + $total3 + $total4 + $total5 + $total6 + $total7 + $total8 + $total9 + $total10;  
        
        
        
        ### Determine number of cannibal and dispersal males and females in all Habitats

        # Find High Dispersal Cannibals
        $lenF1Hc = countHighCan(@femalesA1); $lenM1Hc = countHighCan(@malesA1); 
        $lenF2Hc = countHighCan(@femalesA2); $lenM2Hc = countHighCan(@malesA2); 
        $lenF3Hc = countHighCan(@femalesA3); $lenM3Hc = countHighCan(@malesA3); 
        $lenF4Hc = countHighCan(@femalesA4); $lenM4Hc = countHighCan(@malesA4); 
        $lenF5Hc = countHighCan(@femalesA5); $lenM5Hc = countHighCan(@malesA5); 
        $lenF6Hc = countHighCan(@femalesA6); $lenM6Hc = countHighCan(@malesA6); 
        $lenF7Hc = countHighCan(@femalesA7); $lenM7Hc = countHighCan(@malesA7); 
        $lenF8Hc = countHighCan(@femalesA8); $lenM8Hc = countHighCan(@malesA8); 
        $lenF9Hc = countHighCan(@femalesA9); $lenM9Hc = countHighCan(@malesA9); 
        $lenF10Hc = countHighCan(@femalesA10); $lenM10Hc = countHighCan(@malesA10); 

        # Find Medium Dispersal Cannibals
        $lenF1Mc = countMidCan(@femalesA1); $lenM1Mc = countMidCan(@malesA1); 
        $lenF2Mc = countMidCan(@femalesA2); $lenM2Mc = countMidCan(@malesA2); 
        $lenF3Mc = countMidCan(@femalesA3); $lenM3Mc = countMidCan(@malesA3); 
        $lenF4Mc = countMidCan(@femalesA4); $lenM4Mc = countMidCan(@malesA4); 
        $lenF5Mc = countMidCan(@femalesA5); $lenM5Mc = countMidCan(@malesA5); 
        $lenF6Mc = countMidCan(@femalesA6); $lenM6Mc = countMidCan(@malesA6); 
        $lenF7Mc = countMidCan(@femalesA7); $lenM7Mc = countMidCan(@malesA7); 
        $lenF8Mc = countMidCan(@femalesA8); $lenM8Mc = countMidCan(@malesA8); 
        $lenF9Mc = countMidCan(@femalesA9); $lenM9Mc = countMidCan(@malesA9); 
        $lenF10Mc = countMidCan(@femalesA10); $lenM10Mc = countMidCan(@malesA10); 
        
        # Find Low Dispersal Cannibals
        $lenF1Lc = countLowCan(@femalesA1); $lenM1Lc = countLowCan(@malesA1); 
        $lenF2Lc = countLowCan(@femalesA2); $lenM2Lc = countLowCan(@malesA2); 
        $lenF3Lc = countLowCan(@femalesA3); $lenM3Lc = countLowCan(@malesA3); 
        $lenF4Lc = countLowCan(@femalesA4); $lenM4Lc = countLowCan(@malesA4); 
        $lenF5Lc = countLowCan(@femalesA5); $lenM5Lc = countLowCan(@malesA5); 
        $lenF6Lc = countLowCan(@femalesA6); $lenM6Lc = countLowCan(@malesA6); 
        $lenF7Lc = countLowCan(@femalesA7); $lenM7Lc = countLowCan(@malesA7); 
        $lenF8Lc = countLowCan(@femalesA8); $lenM8Lc = countLowCan(@malesA8); 
        $lenF9Lc = countLowCan(@femalesA9); $lenM9Lc = countLowCan(@malesA9); 
        $lenF10Lc = countLowCan(@femalesA10); $lenM10Lc = countLowCan(@malesA10); 

		$totalA1c = $lenF1Hc + $lenF1Mc + $lenF1Lc + $lenM1Hc + $lenM1Mc + $lenM1Lc;
		$totalA2c = $lenF2Hc + $lenF2Mc + $lenF2Lc + $lenM2Hc + $lenM2Mc + $lenM2Lc;
		$totalA3c = $lenF3Hc + $lenF3Mc + $lenF3Lc + $lenM3Hc + $lenM3Mc + $lenM3Lc;		
        $totalA4c = $lenF4Hc + $lenF4Mc + $lenF4Lc + $lenM4Hc + $lenM4Mc + $lenM4Lc;		
        $totalA5c = $lenF5Hc + $lenF5Mc + $lenF5Lc + $lenM5Hc + $lenM5Mc + $lenM5Lc;		
        $totalA6c = $lenF6Hc + $lenF6Mc + $lenF6Lc + $lenM6Hc + $lenM6Mc + $lenM6Lc;		
        $totalA7c = $lenF7Hc + $lenF7Mc + $lenF7Lc + $lenM7Hc + $lenM7Mc + $lenM7Lc;		
        $totalA8c = $lenF8Hc + $lenF8Mc + $lenF8Lc + $lenM8Hc + $lenM8Mc + $lenM8Lc;		
        $totalA9c = $lenF9Hc + $lenF9Mc + $lenF9Lc + $lenM9Hc + $lenM9Mc + $lenM9Lc;		
        $totalA10c = $lenF10Hc + $lenF10Mc + $lenF10Lc + $lenM10Hc + $lenM10Mc + $lenM10Lc;		
        		
        # Find High Dispersal Vegans
        $lenF1Hv = countHighVegan(@femalesA1); $lenM1Hv = countHighVegan(@malesA1); 
        $lenF2Hv = countHighVegan(@femalesA2); $lenM2Hv = countHighVegan(@malesA2); 
        $lenF3Hv = countHighVegan(@femalesA3); $lenM3Hv = countHighVegan(@malesA3); 
        $lenF4Hv = countHighVegan(@femalesA4); $lenM4Hv = countHighVegan(@malesA4); 
        $lenF5Hv = countHighVegan(@femalesA5); $lenM5Hv = countHighVegan(@malesA5); 
        $lenF6Hv = countHighVegan(@femalesA6); $lenM6Hv = countHighVegan(@malesA6); 
        $lenF7Hv = countHighVegan(@femalesA7); $lenM7Hv = countHighVegan(@malesA7); 
        $lenF8Hv = countHighVegan(@femalesA8); $lenM8Hv = countHighVegan(@malesA8); 
        $lenF9Hv = countHighVegan(@femalesA9); $lenM9Hv = countHighVegan(@malesA9); 
        $lenF10Hv = countHighVegan(@femalesA10); $lenM10Hv = countHighVegan(@malesA10); 

        # Find Medium Dispersal Vegans
        $lenF1Mv = countMidVegan(@femalesA1); $lenM1Mv = countMidVegan(@malesA1); 
        $lenF2Mv = countMidVegan(@femalesA2); $lenM2Mv = countMidVegan(@malesA2); 
        $lenF3Mv = countMidVegan(@femalesA3); $lenM3Mv = countMidVegan(@malesA3); 
        $lenF4Mv = countMidVegan(@femalesA4); $lenM4Mv = countMidVegan(@malesA4); 
        $lenF5Mv = countMidVegan(@femalesA5); $lenM5Mv = countMidVegan(@malesA5); 
        $lenF6Mv = countMidVegan(@femalesA6); $lenM6Mv = countMidVegan(@malesA6); 
        $lenF7Mv = countMidVegan(@femalesA7); $lenM7Mv = countMidVegan(@malesA7); 
        $lenF8Mv = countMidVegan(@femalesA8); $lenM8Mv = countMidVegan(@malesA8); 
        $lenF9Mv = countMidVegan(@femalesA9); $lenM9Mv = countMidVegan(@malesA9); 
        $lenF10Mv = countMidVegan(@femalesA10); $lenM10Mv = countMidVegan(@malesA10); 
        
        
        # Find Low Dispersal Vegans
        $lenF1Lv = countLowVegan(@femalesA1); $lenM1Lv = countLowVegan(@malesA1); 
        $lenF2Lv = countLowVegan(@femalesA2); $lenM2Lv = countLowVegan(@malesA2); 
        $lenF3Lv = countLowVegan(@femalesA3); $lenM3Lv = countLowVegan(@malesA3); 
        $lenF4Lv = countLowVegan(@femalesA4); $lenM4Lv = countLowVegan(@malesA4); 
        $lenF5Lv = countLowVegan(@femalesA5); $lenM5Lv = countLowVegan(@malesA5); 
        $lenF6Lv = countLowVegan(@femalesA6); $lenM6Lv = countLowVegan(@malesA6); 
        $lenF7Lv = countLowVegan(@femalesA7); $lenM7Lv = countLowVegan(@malesA7); 
        $lenF8Lv = countLowVegan(@femalesA8); $lenM8Lv = countLowVegan(@malesA8); 
        $lenF9Lv = countLowVegan(@femalesA9); $lenM9Lv = countLowVegan(@malesA9); 
        $lenF10Lv = countLowVegan(@femalesA10); $lenM10Lv = countLowVegan(@malesA10); 





        # # Female and Male Totals
        ### Totals Female/Male Cannibal and Non-Cannibals for High Dispersal
        $lenF1H = $lenF1Hc + $lenF1Hv; $lenM1H = $lenM1Hc + $lenM1Hv;
        $lenF2H = $lenF2Hc + $lenF2Hv; $lenM2H = $lenM2Hc + $lenM2Hv;
        $lenF3H = $lenF3Hc + $lenF3Hv; $lenM3H = $lenM3Hc + $lenM3Hv;
        $lenF4H = $lenF4Hc + $lenF4Hv; $lenM4H = $lenM4Hc + $lenM4Hv;
        $lenF5H = $lenF5Hc + $lenF5Hv; $lenM5H = $lenM5Hc + $lenM5Hv;
        $lenF6H = $lenF6Hc + $lenF6Hv; $lenM6H = $lenM6Hc + $lenM6Hv;
        $lenF7H = $lenF7Hc + $lenF7Hv; $lenM7H = $lenM7Hc + $lenM7Hv;
        $lenF8H = $lenF8Hc + $lenF8Hv; $lenM8H = $lenM8Hc + $lenM8Hv;
        $lenF9H = $lenF9Hc + $lenF9Hv; $lenM9H = $lenM9Hc + $lenM9Hv;
        $lenF10H = $lenF10Hc + $lenF10Hv; $lenM10H = $lenM10Hc + $lenM10Hv;


        ### Totals Female/Male Cannibal and Non-Cannibals for Medium Dispersal
        $lenF1M = $lenF1Mc + $lenF1Mv; $lenM1M = $lenM1Mc + $lenM1Mv;
        $lenF2M = $lenF2Mc + $lenF2Mv; $lenM2M = $lenM2Mc + $lenM2Mv;
        $lenF3M = $lenF3Mc + $lenF3Mv; $lenM3M = $lenM3Mc + $lenM3Mv;
        $lenF4M = $lenF4Mc + $lenF4Mv; $lenM4M = $lenM4Mc + $lenM4Mv;
        $lenF5M = $lenF5Mc + $lenF5Mv; $lenM5M = $lenM5Mc + $lenM5Mv;
        $lenF6M = $lenF6Mc + $lenF6Mv; $lenM6M = $lenM6Mc + $lenM6Mv;
        $lenF7M = $lenF7Mc + $lenF7Mv; $lenM7M = $lenM7Mc + $lenM7Mv;
        $lenF8M = $lenF8Mc + $lenF8Mv; $lenM8M = $lenM8Mc + $lenM8Mv;
        $lenF9M = $lenF9Mc + $lenF9Mv; $lenM9M = $lenM9Mc + $lenM9Mv;
        $lenF10M = $lenF10Mc + $lenF10Mv; $lenM10M = $lenM10Mc + $lenM10Mv;


        ### Totals Female/Male Cannibal and Non-Cannibals for Low Dispersal
        $lenF1L = $lenF1Lc + $lenF1Lv; $lenM1L = $lenM1Lc + $lenM1Lv;
        $lenF2L = $lenF2Lc + $lenF2Lv; $lenM2L = $lenM2Lc + $lenM2Lv;
        $lenF3L = $lenF3Lc + $lenF3Lv; $lenM3L = $lenM3Lc + $lenM3Lv;
        $lenF4L = $lenF4Lc + $lenF4Lv; $lenM4L = $lenM4Lc + $lenM4Lv;
        $lenF5L = $lenF5Lc + $lenF5Lv; $lenM5L = $lenM5Lc + $lenM5Lv;
        $lenF6L = $lenF6Lc + $lenF6Lv; $lenM6L = $lenM6Lc + $lenM6Lv;
        $lenF7L = $lenF7Lc + $lenF7Lv; $lenM7L = $lenM7Lc + $lenM7Lv;
        $lenF8L = $lenF8Lc + $lenF8Lv; $lenM8L = $lenM8Lc + $lenM8Lv;
        $lenF9L = $lenF9Lc + $lenF9Lv; $lenM9L = $lenM9Lc + $lenM9Lv;
        $lenF10L = $lenF10Lc + $lenF10Lv; $lenM10L = $lenM10Lc + $lenM10Lv;






        # # # Totals for Each Habitat, by Dispersal and Cannibal
        ### Total Cannibals and Non-cannibals for High Dispersal
        $total1Hc = $lenF1Hc + $lenM1Hc; $total1Hv = $lenF1Hv + $lenM1Hv;
        $total2Hc = $lenF2Hc + $lenM2Hc; $total2Hv = $lenF2Hv + $lenM2Hv;
        $total3Hc = $lenF3Hc + $lenM3Hc; $total3Hv = $lenF3Hv + $lenM3Hv;
        $total4Hc = $lenF4Hc + $lenM4Hc; $total4Hv = $lenF4Hv + $lenM4Hv;
        $total5Hc = $lenF5Hc + $lenM5Hc; $total5Hv = $lenF5Hv + $lenM5Hv;
        $total6Hc = $lenF6Hc + $lenM6Hc; $total6Hv = $lenF6Hv + $lenM6Hv;
        $total7Hc = $lenF7Hc + $lenM7Hc; $total7Hv = $lenF7Hv + $lenM7Hv;
        $total8Hc = $lenF8Hc + $lenM8Hc; $total8Hv = $lenF8Hv + $lenM8Hv;
        $total9Hc = $lenF9Hc + $lenM9Hc; $total9Hv = $lenF9Hv + $lenM9Hv;
        $total10Hc = $lenF10Hc + $lenM10Hc; $total10Hv = $lenF10Hv + $lenM10Hv;

        ### Total Cannibals and Non-cannibals for Medium Dispersal
        $total1Mc = $lenF1Mc + $lenM1Mc; $total1Mv = $lenF1Mv + $lenM1Mv;
        $total2Mc = $lenF2Mc + $lenM2Mc; $total2Mv = $lenF2Mv + $lenM2Mv;
        $total3Mc = $lenF3Mc + $lenM3Mc; $total3Mv = $lenF3Mv + $lenM3Mv;
        $total4Mc = $lenF4Mc + $lenM4Mc; $total4Mv = $lenF4Mv + $lenM4Mv;
        $total5Mc = $lenF5Mc + $lenM5Mc; $total5Mv = $lenF5Mv + $lenM5Mv;
        $total6Mc = $lenF6Mc + $lenM6Mc; $total6Mv = $lenF6Mv + $lenM6Mv;
        $total7Mc = $lenF7Mc + $lenM7Mc; $total7Mv = $lenF7Mv + $lenM7Mv;
        $total8Mc = $lenF8Mc + $lenM8Mc; $total8Mv = $lenF8Mv + $lenM8Mv;
        $total9Mc = $lenF9Mc + $lenM9Mc; $total9Mv = $lenF9Mv + $lenM9Mv;
        $total10Mc = $lenF10Mc + $lenM10Mc; $total10Mv = $lenF10Mv + $lenM10Mv;

        ### Total Cannibals and Non-cannibals for Low Dispersal
        $total1Lc = $lenF1Lc + $lenM1Lc; $total1Lv = $lenF1Lv + $lenM1Lv;
        $total2Lc = $lenF2Lc + $lenM2Lc; $total2Lv = $lenF2Lv + $lenM2Lv;
        $total3Lc = $lenF3Lc + $lenM3Lc; $total3Lv = $lenF3Lv + $lenM3Lv;
        $total4Lc = $lenF4Lc + $lenM4Lc; $total4Lv = $lenF4Lv + $lenM4Lv;
        $total5Lc = $lenF5Lc + $lenM5Lc; $total5Lv = $lenF5Lv + $lenM5Lv;
        $total6Lc = $lenF6Lc + $lenM6Lc; $total6Lv = $lenF6Lv + $lenM6Lv;
        $total7Lc = $lenF7Lc + $lenM7Lc; $total7Lv = $lenF7Lv + $lenM7Lv;
        $total8Lc = $lenF8Lc + $lenM8Lc; $total8Lv = $lenF8Lv + $lenM8Lv;
        $total9Lc = $lenF9Lc + $lenM9Lc; $total9Lv = $lenF9Lv + $lenM9Lv;
        $total10Lc = $lenF10Lc + $lenM10Lc; $total10Lv = $lenF10Lv + $lenM10Lv;
 
        
        ##################################################
        ########## Determine allele frequencies ##########
        ##################################################
        
        

           	
    	### Calculate proportion of Species A in each habitat 
    	$propA1 = calcSpeciesProportion($totalA1, $total1); 
    	$propA2 = calcSpeciesProportion($totalA2, $total2); 
    	$propA3 = calcSpeciesProportion($totalA3, $total3); 
    	$propA4 = calcSpeciesProportion($totalA4, $total4); 
    	$propA5 = calcSpeciesProportion($totalA5, $total5); 
    	$propA6 = calcSpeciesProportion($totalA6, $total6); 
    	$propA7 = calcSpeciesProportion($totalA7, $total7); 
    	$propA8 = calcSpeciesProportion($totalA8, $total8); 
    	$propA9 = calcSpeciesProportion($totalA9, $total9); 
    	$propA10 = calcSpeciesProportion($totalA10, $total10); 
        
        #total proportion of Species A in overall simulation (all 9 habitats)
        $propATotal = calcSpeciesProportion($totalA, $total); 
        
    	### Calculate proportion of Species A in each habitat 
    	$propB1 = calcSpeciesProportion($totalB1, $total1); 
    	$propB2 = calcSpeciesProportion($totalB2, $total2); 
    	$propB3 = calcSpeciesProportion($totalB3, $total3); 
    	$propB4 = calcSpeciesProportion($totalB4, $total4); 
    	$propB5 = calcSpeciesProportion($totalB5, $total5); 
    	$propB6 = calcSpeciesProportion($totalB6, $total6); 
    	$propB7 = calcSpeciesProportion($totalB7, $total7); 
    	$propB8 = calcSpeciesProportion($totalB8, $total8); 
    	$propB9 = calcSpeciesProportion($totalB9, $total9); 
    	$propB10 = calcSpeciesProportion($totalB10, $total10); 
        
        #total proportion of Species B in overall simulation (all 9 habitats)
        $propBTotal = calcSpeciesProportion($totalB, $total); 
        
        
        # Calculate proportion of cannibals within each species in each habitat
        $propA1c = calcSpeciesProportion($totalA1c,$totalA1); 
        $propA2c = calcSpeciesProportion($totalA2c,$totalA2); 
        $propA3c = calcSpeciesProportion($totalA3c,$totalA3); 
        $propA4c = calcSpeciesProportion($totalA4c,$totalA4); 
        $propA5c = calcSpeciesProportion($totalA5c,$totalA5); 
        $propA6c = calcSpeciesProportion($totalA6c,$totalA6); 
        $propA7c = calcSpeciesProportion($totalA7c,$totalA7); 
        $propA8c = calcSpeciesProportion($totalA8c,$totalA8); 
        $propA9c = calcSpeciesProportion($totalA9c,$totalA9); 
        $propA10c = calcSpeciesProportion($totalA10c,$totalA10); 
        
        $propAcTotal = calcSpeciesProportion($totalAc,$totalA); 
        
        $propB1c = calcSpeciesProportion($totalB1c,$totalB1); 
        $propB2c = calcSpeciesProportion($totalB2c,$totalB2); 
        $propB3c = calcSpeciesProportion($totalB3c,$totalB3); 
        $propB4c = calcSpeciesProportion($totalB4c,$totalB4); 
        $propB5c = calcSpeciesProportion($totalB5c,$totalB5); 
        $propB6c = calcSpeciesProportion($totalB6c,$totalB6); 
        $propB7c = calcSpeciesProportion($totalB7c,$totalB7); 
        $propB8c = calcSpeciesProportion($totalB8c,$totalB8); 
        $propB9c = calcSpeciesProportion($totalB9c,$totalB9); 
        $propB10c = calcSpeciesProportion($totalB10c,$totalB10); 
        
        $propBcTotal = calcSpeciesProportion($totalBc,$totalB); 
               
        
        # Calculate proportion of cannibal females within Species A and B in each habitat 

        
        # Total proportion of cannibal females in each Species in all habitats 
        
    	
       
        ### Increase generation by 1
        $gen++;
#print "propa1: $lenF1Hc\n";        
        # Rounded values of proportions of each Species 
        $rpropA1 = sprintf "%.3f", $propA1; 
        $rpropA2 = sprintf "%.3f", $propA2; 
        $rpropA3 = sprintf "%.3f", $propA3; 
        $rpropA4 = sprintf "%.3f", $propA4; 
        $rpropA5 = sprintf "%.3f", $propA5; 
        $rpropA6 = sprintf "%.3f", $propA6; 
        $rpropA7 = sprintf "%.3f", $propA7; 
        $rpropA8 = sprintf "%.3f", $propA8; 
        $rpropA9 = sprintf "%.3f", $propA9; 
        $rpropA10 = sprintf "%.3f", $propA10; 
        $rpropATotal = sprintf "%.3f", $propATotal; 
        
        $rpropB1 = sprintf "%.3f", $propB1; 
        $rpropB2 = sprintf "%.3f", $propB2; 
        $rpropB3 = sprintf "%.3f", $propB3; 
        $rpropB4 = sprintf "%.3f", $propB4; 
        $rpropB5 = sprintf "%.3f", $propB5; 
        $rpropB6 = sprintf "%.3f", $propB6; 
        $rpropB7 = sprintf "%.3f", $propB7; 
        $rpropB8 = sprintf "%.3f", $propB8; 
        $rpropB9 = sprintf "%.3f", $propB9; 
        $rpropB10 = sprintf "%.3f", $propB10; 
        $rpropBTotal = sprintf "%.3f", $propBTotal; 
        
        $rpropA1c = sprintf "%.3f", $propA1c; 
        $rpropA2c = sprintf "%.3f", $propA2c; 
        $rpropA3c = sprintf "%.3f", $propA3c; 
        $rpropA4c = sprintf "%.3f", $propA4c; 
        $rpropA5c = sprintf "%.3f", $propA5c; 
        $rpropA6c = sprintf "%.3f", $propA6c; 
        $rpropA7c = sprintf "%.3f", $propA7c; 
        $rpropA8c = sprintf "%.3f", $propA8c; 
        $rpropA9c = sprintf "%.3f", $propA9c; 
        $rpropA10c = sprintf "%.3f", $propA10c; 
        $rpropAcTotal = sprintf "%.3f", $propAcTotal; 

        $rpropB1c = sprintf "%.3f", $propB1c; 
        $rpropB2c = sprintf "%.3f", $propB2c; 
        $rpropB3c = sprintf "%.3f", $propB3c; 
        $rpropB4c = sprintf "%.3f", $propB4c; 
        $rpropB5c = sprintf "%.3f", $propB5c; 
        $rpropB6c = sprintf "%.3f", $propB6c; 
        $rpropB7c = sprintf "%.3f", $propB7c; 
        $rpropB8c = sprintf "%.3f", $propB8c; 
        $rpropB9c = sprintf "%.3f", $propB9c; 
        $rpropB10c = sprintf "%.3f", $propB10c; 
        $rpropBcTotal = sprintf "%.3f", $propBcTotal; 
        
               
        ### Terminal Output ###  
        
        #Population information for Habitats 1-5
        #print "Simulation $simno\tGeneration $gen\n"; 
        # print "Mutation $mutrt; Migration A $migrtA; Migration B $migrtB; Cancost $cancost_rand; Can_offspring $num_cannibal_offspring; Capacity $carrying_capacity\n"; 
        # print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n"; 
        # print "Habitat 1\t\t\tHabitat 2\t\t\tHabitat 3\t\tHabitat 4\t\tHabitat 5\n";
        # 
        # print "Species A\t\t\tSpecies A\t\t\tSpecies A\t\tSpecies A\t\tSpecies A\n"; 
        # print "total female : $lenFA1\t\ttotal female : $lenFA2\t\ttotal female : $lenFA3\ttotal female : $lenFA4\ttotal female : $lenFA5\n";
        # print "total male   : $lenMA1\t\ttotal male   : $lenMA2\t\ttotal male   : $lenMA3\ttotal male : $lenMA4\ttotal male : $lenMA5\n";
        # print "total : $totalA1\t\t\ttotal : $totalA2\t\t\ttotal : $totalA3\t\ttotal : $totalA4\t\ttotal : $totalA5\n\n";
        # 
        # print "Species B\t\t\tSpecies B\t\t\tSpecies B\t\tSpecies B\t\tSpecies B\n"; 
        # print "total female : $lenFB1\t\ttotal female : $lenFB2\t\ttotal female : $lenFB3\ttotal female : $lenFB4\ttotal female : $lenFB5\n";
        # print "total male   : $lenMB1\t\ttotal male   : $lenMB2\t\ttotal male   : $lenMB3\ttotal male : $lenMB4\ttotal male : $lenMB5\n";
        # print "total : $totalB1\t\t\ttotal : $totalB2\t\t\ttotal : $totalB3\t\ttotal : $totalB4\t\ttotal : $totalB5\n\n";
        # print "===================================================================================================================\n\n";
        # 
        # 
        # #Population information for Habitats 6-10
        # print "Habitat 6\t\t\tHabitat 7\t\t\tHabitat 8\t\tHabitat 9\t\tHabitat 10\n";
        # 
        # print "Species A\t\t\tSpecies A\t\t\tSpecies A\t\tSpecies A\t\tSpecies A\n"; 
        # print "total female : $lenFA6\t\ttotal female : $lenFA7\t\ttotal female : $lenFA8\ttotal female : $lenFA9\ttotal female : $lenFA10\n";
        # print "total male   : $lenMA6\t\ttotal male   : $lenMA7\t\ttotal male   : $lenMA8\ttotal male : $lenMA9\ttotal male : $lenMA10\n";
        # print "total : $totalA6\t\t\ttotal : $totalA7\t\t\ttotal : $totalA8\t\ttotal : $totalA9\t\ttotal : $totalA10\n\n";
        # 
        # print "Species B\t\t\tSpecies B\t\t\tSpecies B\t\tSpecies B\t\tSpecies B\n"; 
        # print "total female : $lenFB6\t\ttotal female : $lenFB7\t\ttotal female : $lenFB8\ttotal female : $lenFB9\ttotal female : $lenFB10\n";
        # print "total male   : $lenMB6\t\ttotal male   : $lenMB7\t\ttotal male   : $lenMB8\ttotal male : $lenMB9\ttotal male : $lenMB10\n";
        # print "total : $totalB6\t\t\ttotal : $totalB7\t\t\ttotal : $totalB8\t\ttotal : $totalB9\t\ttotal : $totalB10\n\n";
        # print "===================================================================================================================\n";
        # 
        # 
       # 
       
       
        ##################################################
        ################ Create offspring ################
        ##################################################
        
        ### Create arrays to hold offspring
        # Name key:
            # Format: (F/M)offspring(A/B)(#)
            # F/M = fe/male offspring
            # A/B = Species A/Species B
            # # = habitat number
        
        # Female and male offspring of Species A in Habitats 1-9
     
        @FoffspringA1 = (); @MoffspringA1 = (); 
        @FoffspringA2 = (); @MoffspringA2 = (); 
        @FoffspringA3 = (); @MoffspringA3 = (); 
        @FoffspringA4 = (); @MoffspringA4 = (); 
        @FoffspringA5 = (); @MoffspringA5 = (); 
        @FoffspringA6 = (); @MoffspringA6 = (); 
        @FoffspringA7 = (); @MoffspringA7 = (); 
        @FoffspringA8 = (); @MoffspringA8 = (); 
        @FoffspringA9 = (); @MoffspringA9 = (); 
        @FoffspringA10 = (); @MoffspringA10 = (); 
        
        # Female and male offspring of Species B in Habitats 1-9
     
        @FoffspringB1 = (); @MoffspringB1 = (); 
        @FoffspringB2 = (); @MoffspringB2 = (); 
        @FoffspringB3 = (); @MoffspringB3 = (); 
        @FoffspringB4 = (); @MoffspringB4 = (); 
        @FoffspringB5 = (); @MoffspringB5 = (); 
        @FoffspringB6 = (); @MoffspringB6 = (); 
        @FoffspringB7 = (); @MoffspringB7 = (); 
        @FoffspringB8 = (); @MoffspringB8 = (); 
        @FoffspringB9 = (); @MoffspringB9 = (); 
        @FoffspringB10 = (); @MoffspringB10 = (); 
        
               
        ### Calculate maximum offspring that the habitat (both speciess) can produce that generation--will always be equal to carrying capacity  

        $max_offspring1 = calcCarryingCapacity(); #$lenFA1, $lenMA1, $lenFB1, $lenMB1); #, $k1);
        $max_offspring2 = calcCarryingCapacity(); #$lenFA2, $lenMA2, $lenFB2, $lenMB2); #, $k2);
        $max_offspring3 = calcCarryingCapacity(); #$lenFA3, $lenMA3, $lenFB3, $lenMB3); #, $k3);
        $max_offspring4 = calcCarryingCapacity(); #$lenFA4, $lenMA4, $lenFB4, $lenMB4); #, $k4);
        $max_offspring5 = calcCarryingCapacity(); #$lenFA5, $lenMA5, $lenFB5, $lenMB5); #, $k5);
        $max_offspring6 = calcCarryingCapacity(); #$lenFA6, $lenMA6, $lenFB6, $lenMB6); #, $k6);
        $max_offspring7 = calcCarryingCapacity(); #$lenFA7, $lenMA7, $lenFB7, $lenMB7); #, $k7);
        $max_offspring8 = calcCarryingCapacity(); #$lenFA8, $lenMA8, $lenFB8, $lenMB8); #, $k8);
        $max_offspring9 = calcCarryingCapacity(); #$lenFA9, $lenMA9, $lenFB9, $lenMB9); #, $k9);
        $max_offspring10 = calcCarryingCapacity(); #$lenFA10, $lenMA10, $lenFB10, $lenMB10); #, $k10);
        

        
        ########################################
        ###### Offspring subroutine calls ######
        ########################################
        
        # Produce the next generation for both species at each habitat.   
        
        @female_trackerA1 = (); @male_trackerA1 = (); 
        @female_trackerA2 = (); @male_trackerA2 = (); 
        @female_trackerA3 = (); @male_trackerA3 = (); 
        @female_trackerA4 = (); @male_trackerA4 = (); 
        @female_trackerA5 = (); @male_trackerA5 = (); 
        @female_trackerA6 = (); @male_trackerA6 = (); 
        @female_trackerA7 = (); @male_trackerA7 = (); 
        @female_trackerA8 = (); @male_trackerA8 = (); 
        @female_trackerA9 = (); @male_trackerA9 = (); 
        @female_trackerA10 = (); @male_trackerA10 = (); 

        @female_trackerB1 = (); @male_trackerB1 = (); 
        @female_trackerB2 = (); @male_trackerB2 = (); 
        @female_trackerB3 = (); @male_trackerB3 = (); 
        @female_trackerB4 = (); @male_trackerB4 = (); 
        @female_trackerB5 = (); @male_trackerB5 = (); 
        @female_trackerB6 = (); @male_trackerB6 = (); 
        @female_trackerB7 = (); @male_trackerB7 = (); 
        @female_trackerB8 = (); @male_trackerB8 = (); 
        @female_trackerB9 = (); @male_trackerB9 = (); 
        @female_trackerB10 = (); @male_trackerB10 = (); 
        
        makeNextGenCompetingSpecies($lenF1H, $lenM1H, $lenF1M, $lenM1M, $lenF1L, $lenM1L, $max_offspring1, @femalesA1, @malesA1, @femalesB1, @malesB1, @FoffspringA1, @MoffspringA1, @FoffspringB1, @MoffspringB1, @female_trackerA1, @male_trackerA1, @female_trackerB1, @male_trackerB1); 
        makeNextGenCompetingSpecies($lenF2H, $lenM2H, $lenF2M, $lenM2M, $lenF2L, $lenM2L, $max_offspring2, @femalesA2, @malesA2, @femalesB2, @malesB2, @FoffspringA2, @MoffspringA2, @FoffspringB2, @MoffspringB2, @female_trackerA2, @male_trackerA2, @female_trackerB2, @male_trackerB2); 
        makeNextGenCompetingSpecies($lenF3H, $lenM3H, $lenF3M, $lenM3M, $lenF3L, $lenM3L, $max_offspring3, @femalesA3, @malesA3, @femalesB3, @malesB3, @FoffspringA3, @MoffspringA3, @FoffspringB3, @MoffspringB3, @female_trackerA3, @male_trackerA3, @female_trackerB3, @male_trackerB3); 
        makeNextGenCompetingSpecies($lenF4H, $lenM4H, $lenF4M, $lenM4M, $lenF4L, $lenM4L, $max_offspring4, @femalesA4, @malesA4, @femalesB4, @malesB4, @FoffspringA4, @MoffspringA4, @FoffspringB4, @MoffspringB4, @female_trackerA4, @male_trackerA4, @female_trackerB4, @male_trackerB4); 
        makeNextGenCompetingSpecies($lenF5H, $lenM5H, $lenF5M, $lenM5M, $lenF5L, $lenM5L, $max_offspring5, @femalesA5, @malesA5, @femalesB5, @malesB5, @FoffspringA5, @MoffspringA5, @FoffspringB5, @MoffspringB5, @female_trackerA5, @male_trackerA5, @female_trackerB5, @male_trackerB5); 
        makeNextGenCompetingSpecies($lenF6H, $lenM6H, $lenF6M, $lenM6M, $lenF6L, $lenM6L, $max_offspring6, @femalesA6, @malesA6, @femalesB6, @malesB6, @FoffspringA6, @MoffspringA6, @FoffspringB6, @MoffspringB6, @female_trackerA6, @male_trackerA6, @female_trackerB6, @male_trackerB6); 
        makeNextGenCompetingSpecies($lenF7H, $lenM7H, $lenF7M, $lenM7M, $lenF7L, $lenM7L, $max_offspring7, @femalesA7, @malesA7, @femalesB7, @malesB7, @FoffspringA7, @MoffspringA7, @FoffspringB7, @MoffspringB7, @female_trackerA7, @male_trackerA7, @female_trackerB7, @male_trackerB7); 
        makeNextGenCompetingSpecies($lenF8H, $lenM8H, $lenF8M, $lenM8M, $lenF8L, $lenM8L, $max_offspring8, @femalesA8, @malesA8, @femalesB8, @malesB8, @FoffspringA8, @MoffspringA8, @FoffspringB8, @MoffspringB8, @female_trackerA8, @male_trackerA8, @female_trackerB8, @male_trackerB8); 
        makeNextGenCompetingSpecies($lenF9H, $lenM9H, $lenF9M, $lenM9M, $lenF9L, $lenM9L, $max_offspring9, @femalesA9, @malesA9, @femalesB9, @malesB9, @FoffspringA9, @MoffspringA9, @FoffspringB9, @MoffspringB9, @female_trackerA9, @male_trackerA9, @female_trackerB9, @male_trackerB9); 
        makeNextGenCompetingSpecies($lenF10H, $lenM10H, $lenF10M, $lenM10M, $lenF10L, $lenM10L, $max_offspring10, @femalesA10, @malesA10, @femalesB10, @malesB10, @FoffspringA10, @MoffspringA10, @FoffspringB10, @MoffspringB10, @female_trackerA10, @male_trackerA10, @female_trackerB10, @male_trackerB10); 
               
        # print "Females: "; 
        # foreach(@female_trackerA1){ print "$_\t"; }
        # print "\n"; 

        ##################################################
        #################### Migration ###################
        ##################################################
        
  
            ### Create array to store new locations of individuals after migration
            # Name key:
                # Format: (F/M)migrants(A/B)(#)
                # F/M = fe/male migrants
                # A/B = Species A/Species B individuals
                # # = habitat number
        
            ### these are the destination arrays; e.g. any Species female that gets put in @FmigrantsA7 will move to habitat 7 and be counted among that female Species A population

            @FmigrantsA1 = (); @MmigrantsA1 = (); 
            @FmigrantsA2 = (); @MmigrantsA2 = (); 
            @FmigrantsA3 = (); @MmigrantsA3 = (); 
            @FmigrantsA4 = (); @MmigrantsA4 = (); 
            @FmigrantsA5 = (); @MmigrantsA5 = (); 
            @FmigrantsA6 = (); @MmigrantsA6 = (); 
            @FmigrantsA7 = (); @MmigrantsA7 = (); 
            @FmigrantsA8 = (); @MmigrantsA8 = (); 
            @FmigrantsA9 = (); @MmigrantsA9 = (); 
            @FmigrantsA10 = (); @MmigrantsA10 = (); 
            
            
            @FmigrantsB1 = (); @MmigrantsB1 = (); 
            @FmigrantsB2 = (); @MmigrantsB2 = (); 
            @FmigrantsB3 = (); @MmigrantsB3 = (); 
            @FmigrantsB4 = (); @MmigrantsB4 = (); 
            @FmigrantsB5 = (); @MmigrantsB5 = (); 
            @FmigrantsB6 = (); @MmigrantsB6 = (); 
            @FmigrantsB7 = (); @MmigrantsB7 = (); 
            @FmigrantsB8 = (); @MmigrantsB8 = (); 
            @FmigrantsB9 = (); @MmigrantsB9 = ();             
            @FmigrantsB10 = (); @MmigrantsB10 = ();             
        
        
            ### Make offspring migrate
            # Habitat 1
            migration(1, @femalesA1, @female_trackerA1, @male_trackerA1, @FmigrantsA1, @FmigrantsA2, @FmigrantsA3, @FmigrantsA4, @FmigrantsA5, @FmigrantsA6, @FmigrantsA7, @FmigrantsA8, @FmigrantsA9, @FmigrantsA10, 
                                @MmigrantsA1, @MmigrantsA2, @MmigrantsA3, @MmigrantsA4, @MmigrantsA5, @MmigrantsA6, @MmigrantsA7, @MmigrantsA8, @MmigrantsA9, @MmigrantsA10, 
                                @FoffspringA1, @MoffspringA1);
                                
            migration(2, @femalesA2, @female_trackerA2, @male_trackerA2, @FmigrantsA1, @FmigrantsA2, @FmigrantsA3, @FmigrantsA4, @FmigrantsA5, @FmigrantsA6, @FmigrantsA7, @FmigrantsA8, @FmigrantsA9, @FmigrantsA10, 
                                @MmigrantsA1, @MmigrantsA2, @MmigrantsA3, @MmigrantsA4, @MmigrantsA5, @MmigrantsA6, @MmigrantsA7, @MmigrantsA8, @MmigrantsA9, @MmigrantsA10, 
                                @FoffspringA2, @MoffspringA2);
                                
            migration(3, @femalesA3, @female_trackerA3, @male_trackerA3, @FmigrantsA1, @FmigrantsA2, @FmigrantsA3, @FmigrantsA4, @FmigrantsA5, @FmigrantsA6, @FmigrantsA7, @FmigrantsA8, @FmigrantsA9, @FmigrantsA10, 
                                @MmigrantsA1, @MmigrantsA2, @MmigrantsA3, @MmigrantsA4, @MmigrantsA5, @MmigrantsA6, @MmigrantsA7, @MmigrantsA8, @MmigrantsA9, @MmigrantsA10, 
                                @FoffspringA3, @MoffspringA3);
                                
            migration(4, @femalesA4, @female_trackerA4, @male_trackerA4, @FmigrantsA1, @FmigrantsA2, @FmigrantsA3, @FmigrantsA4, @FmigrantsA5, @FmigrantsA6, @FmigrantsA7, @FmigrantsA8, @FmigrantsA9, @FmigrantsA10, 
                                @MmigrantsA1, @MmigrantsA2, @MmigrantsA3, @MmigrantsA4, @MmigrantsA5, @MmigrantsA6, @MmigrantsA7, @MmigrantsA8, @MmigrantsA9, @MmigrantsA10, 
                                @FoffspringA4, @MoffspringA4);
                                
            migration(5, @femalesA5, @female_trackerA5, @male_trackerA5, @FmigrantsA1, @FmigrantsA2, @FmigrantsA3, @FmigrantsA4, @FmigrantsA5, @FmigrantsA6, @FmigrantsA7, @FmigrantsA8, @FmigrantsA9, @FmigrantsA10, 
                                @MmigrantsA1, @MmigrantsA2, @MmigrantsA3, @MmigrantsA4, @MmigrantsA5, @MmigrantsA6, @MmigrantsA7, @MmigrantsA8, @MmigrantsA9, @MmigrantsA10, 
                                @FoffspringA5, @MoffspringA5);
                                
            migration(6, @femalesA6, @female_trackerA6, @male_trackerA6, @FmigrantsA1, @FmigrantsA2, @FmigrantsA3, @FmigrantsA4, @FmigrantsA5, @FmigrantsA6, @FmigrantsA7, @FmigrantsA8, @FmigrantsA9, @FmigrantsA10, 
                                @MmigrantsA1, @MmigrantsA2, @MmigrantsA3, @MmigrantsA4, @MmigrantsA5, @MmigrantsA6, @MmigrantsA7, @MmigrantsA8, @MmigrantsA9, @MmigrantsA10, 
                                @FoffspringA6, @MoffspringA6);
                                
            migration(7, @femalesA7, @female_trackerA7, @male_trackerA7, @FmigrantsA1, @FmigrantsA2, @FmigrantsA3, @FmigrantsA4, @FmigrantsA5, @FmigrantsA6, @FmigrantsA7, @FmigrantsA8, @FmigrantsA9, @FmigrantsA10, 
                                @MmigrantsA1, @MmigrantsA2, @MmigrantsA3, @MmigrantsA4, @MmigrantsA5, @MmigrantsA6, @MmigrantsA7, @MmigrantsA8, @MmigrantsA9, @MmigrantsA10, 
                                @FoffspringA7, @MoffspringA7);
                                
            migration(8, @femalesA8, @female_trackerA8, @male_trackerA8, @FmigrantsA1, @FmigrantsA2, @FmigrantsA3, @FmigrantsA4, @FmigrantsA5, @FmigrantsA6, @FmigrantsA7, @FmigrantsA8, @FmigrantsA9, @FmigrantsA10, 
                                @MmigrantsA1, @MmigrantsA2, @MmigrantsA3, @MmigrantsA4, @MmigrantsA5, @MmigrantsA6, @MmigrantsA7, @MmigrantsA8, @MmigrantsA9, @MmigrantsA10, 
                                @FoffspringA8, @MoffspringA8);
                                
            migration(9, @femalesA9, @female_trackerA9, @male_trackerA9, @FmigrantsA1, @FmigrantsA2, @FmigrantsA3, @FmigrantsA4, @FmigrantsA5, @FmigrantsA6, @FmigrantsA7, @FmigrantsA8, @FmigrantsA9, @FmigrantsA10, 
                                @MmigrantsA1, @MmigrantsA2, @MmigrantsA3, @MmigrantsA4, @MmigrantsA5, @MmigrantsA6, @MmigrantsA7, @MmigrantsA8, @MmigrantsA9, @MmigrantsA10, 
                                @FoffspringA9, @MoffspringA9);
                                
            migration(10, @femalesA10, @female_trackerA10, @male_trackerA10, @FmigrantsA1, @FmigrantsA2, @FmigrantsA3, @FmigrantsA4, @FmigrantsA5, @FmigrantsA6, @FmigrantsA7, @FmigrantsA8, @FmigrantsA9, @FmigrantsA10, 
                                @MmigrantsA1, @MmigrantsA2, @MmigrantsA3, @MmigrantsA4, @MmigrantsA5, @MmigrantsA6, @MmigrantsA7, @MmigrantsA8, @MmigrantsA9, @MmigrantsA10, 
                                @FoffspringA10, @MoffspringA10);
                                
                                              
        
        
        ##################################################
        ######### Reset and re-fill adult arrays #########
        ##################################################
        
        @femalesA1 = (); @malesA1 = (); 
        @femalesA2 = (); @malesA2 = (); 
        @femalesA3 = (); @malesA3 = (); 
        @femalesA4 = (); @malesA4 = (); 
        @femalesA5 = (); @malesA5 = (); 
        @femalesA6 = (); @malesA6 = (); 
        @femalesA7 = (); @malesA7 = ();
        @femalesA8 = (); @malesA8 = ();
        @femalesA9 = (); @malesA9 = ();        
        @femalesA10 = (); @malesA10 = ();        
        
        @femalesB1 = (); @malesB1 = (); 
        @femalesB2 = (); @malesB2 = (); 
        @femalesB3 = (); @malesB3 = (); 
        @femalesB4 = (); @malesB4 = (); 
        @femalesB5 = (); @malesB5 = (); 
        @femalesB6 = (); @malesB6 = (); 
        @femalesB7 = (); @malesB7 = ();
        @femalesB8 = (); @malesB8 = ();
        @femalesB9 = (); @malesB9 = ();
        @femalesB10 = (); @malesB10 = ();
        
        ### Fill arrays -- after migration, use migrant arrays; before, use offspring arrays
            shiftPopulation(@femalesA1, @FmigrantsA1, @malesA1, @MmigrantsA1);
            shiftPopulation(@femalesA2, @FmigrantsA2, @malesA2, @MmigrantsA2);
            shiftPopulation(@femalesA3, @FmigrantsA3, @malesA3, @MmigrantsA3);
            shiftPopulation(@femalesA4, @FmigrantsA4, @malesA4, @MmigrantsA4);
            shiftPopulation(@femalesA5, @FmigrantsA5, @malesA5, @MmigrantsA5);
            shiftPopulation(@femalesA6, @FmigrantsA6, @malesA6, @MmigrantsA6);
            shiftPopulation(@femalesA7, @FmigrantsA7, @malesA7, @MmigrantsA7);
            shiftPopulation(@femalesA8, @FmigrantsA8, @malesA8, @MmigrantsA8);
            shiftPopulation(@femalesA9, @FmigrantsA9, @malesA9, @MmigrantsA9);            
            shiftPopulation(@femalesA10, @FmigrantsA10, @malesA10, @MmigrantsA10);            
            
            shiftPopulation(@femalesB1, @FmigrantsB1, @malesB1, @MmigrantsB1);
            shiftPopulation(@femalesB2, @FmigrantsB2, @malesB2, @MmigrantsB2);
            shiftPopulation(@femalesB3, @FmigrantsB3, @malesB3, @MmigrantsB3);
            shiftPopulation(@femalesB4, @FmigrantsB4, @malesB4, @MmigrantsB4);
            shiftPopulation(@femalesB5, @FmigrantsB5, @malesB5, @MmigrantsB5);
            shiftPopulation(@femalesB6, @FmigrantsB6, @malesB6, @MmigrantsB6);
            shiftPopulation(@femalesB7, @FmigrantsB7, @malesB7, @MmigrantsB7);
            shiftPopulation(@femalesB8, @FmigrantsB8, @malesB8, @MmigrantsB8);
            shiftPopulation(@femalesB9, @FmigrantsB9, @malesB9, @MmigrantsB9);
            shiftPopulation(@femalesB10, @FmigrantsB10, @malesB10, @MmigrantsB10);
  
        
        
        #print(FH2 "$totalA1c,$totalA1n,$totalA2c,$totalA2n,$totalA3c,$totalA3n,$totalA4c,$totalA4n,$totalA5c,$totalA5n,$totalA6c,$totalA6n,$totalA7c,$totalA7n,$totalA8c,$totalA8n,$totalA9c,$totalA9n,$totalA10c,$totalA10n,$totalB1c,$totalB1n,$totalB2c,$totalB2n,$totalB3c,$totalB3n,$totalB4c,$totalB4n,$totalB5c,$totalB5n,$totalB6c,$totalB6n,$totalB7c,$totalB7n,$totalB8c,$totalB8n,$totalB9c,$totalB9n,$totalB10c,$totalB10n,$chance_cannibalism_cost,$chance_cannibalism_offspring,$gen,$simno\n"); 
        #print(FH3 "$totalA1c,$totalA1n,$totalA2c,$totalA2n,$totalA3c,$totalA3n,$totalA4c,$totalA4n,$totalA5c,$totalA5n,$totalA6c,$totalA6n,$totalA7c,$totalA7n,$totalA8c,$totalA8n,$totalA9c,$totalA9n,$totalA10c,$totalA10n,$totalB1c,$totalB1n,$totalB2c,$totalB2n,$totalB3c,$totalB3n,$totalB4c,$totalB4n,$totalB5c,$totalB5n,$totalB6c,$totalB6n,$totalB7c,$totalB7n,$totalB8c,$totalB8n,$totalB9c,$totalB9n,$totalB10c,$totalB10n,$gen,$simno\n"); 

	if($gen==$endgen){
        print(FH4 "$gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,1,high,no,$total1Hv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,1,high,yes,$total1Hc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,2,high,no,$total2Hv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,2,high,yes,$total2Hc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,3,high,no,$total3Hv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,3,high,yes,$total3Hc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,4,high,no,$total4Hv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,4,high,yes,$total4Hc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,5,high,no,$total5Hv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,5,high,yes,$total5Hc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,6,high,no,$total6Hv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,6,high,yes,$total6Hc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,7,high,no,$total7Hv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,7,high,yes,$total7Hc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,8,high,no,$total8Hv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,8,high,yes,$total8Hc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,9,high,no,$total9Hv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,9,high,yes,$total9Hc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,10,high,no,$total10Hv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,10,high,yes,$total10Hc,                
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,1,mid,no,$total1Mv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,1,mid,yes,$total1Mc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,2,mid,no,$total2Mv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,2,mid,yes,$total2Mc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,3,mid,no,$total3Mv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,3,mid,yes,$total3Mc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,4,mid,no,$total4Mv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,4,mid,yes,$total4Mc,        
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,5,mid,no,$total5Mv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,5,mid,yes,$total5Mc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,6,mid,no,$total6Mv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,6,mid,yes,$total6Mc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,7,mid,no,$total7Mv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,7,mid,yes,$total7Mc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,8,mid,no,$total8Mv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,8,mid,yes,$total8Mc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,9,mid,no,$total9Mv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,9,mid,yes,$total9Mc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,10,mid,no,$total10Mv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,10,mid,yes,$total10Mc,    
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,1,low,no,$total1Lv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,1,low,yes,$total1Lc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,2,low,no,$total2Lv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,2,low,yes,$total2Lc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,3,low,no,$total3Lv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,3,low,yes,$total3Lc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,4,low,no,$total4Lv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,4,low,yes,$total4Lc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,5,low,no,$total5Lv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,5,low,yes,$total5Lc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,6,low,no,$total6Lv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,6,low,yes,$total6Lc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,7,low,no,$total1Lv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,7,low,yes,$total1Lc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,8,low,no,$total8Lv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,8,low,yes,$total8Lc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,9,low,no,$total9Lv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,9,low,yes,$total9Lc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,10,low,no,$total10Lv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,10,low,yes,$total10Lc\n");
		
print(FH4d "$gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,1,high,no,$total1Hv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,1,high,yes,$total1Hc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,2,high,no,$total2Hv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,2,high,yes,$total2Hc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,3,high,no,$total3Hv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,3,high,yes,$total3Hc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,4,high,no,$total4Hv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,4,high,yes,$total4Hc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,5,high,no,$total5Hv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,5,high,yes,$total5Hc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,6,high,no,$total6Hv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,6,high,yes,$total6Hc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,7,high,no,$total7Hv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,7,high,yes,$total7Hc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,8,high,no,$total8Hv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,8,high,yes,$total8Hc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,9,high,no,$total9Hv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,9,high,yes,$total9Hc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,10,high,no,$total10Hv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,10,high,yes,$total10Hc,                
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,1,mid,no,$total1Mv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,1,mid,yes,$total1Mc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,2,mid,no,$total2Mv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,2,mid,yes,$total2Mc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,3,mid,no,$total3Mv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,3,mid,yes,$total3Mc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,4,mid,no,$total4Mv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,4,mid,yes,$total4Mc,        
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,5,mid,no,$total5Mv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,5,mid,yes,$total5Mc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,6,mid,no,$total6Mv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,6,mid,yes,$total6Mc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,7,mid,no,$total7Mv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,7,mid,yes,$total7Mc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,8,mid,no,$total8Mv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,8,mid,yes,$total8Mc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,9,mid,no,$total9Mv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,9,mid,yes,$total9Mc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,10,mid,no,$total10Mv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,10,mid,yes,$total10Mc,    
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,1,low,no,$total1Lv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,1,low,yes,$total1Lc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,2,low,no,$total2Lv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,2,low,yes,$total2Lc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,3,low,no,$total3Lv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,3,low,yes,$total3Lc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,4,low,no,$total4Lv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,4,low,yes,$total4Lc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,5,low,no,$total5Lv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,5,low,yes,$total5Lc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,6,low,no,$total6Lv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,6,low,yes,$total6Lc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,7,low,no,$total1Lv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,7,low,yes,$total1Lc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,8,low,no,$total8Lv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,8,low,yes,$total8Lc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,9,low,no,$total9Lv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,9,low,yes,$total9Lc,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,10,low,no,$total10Lv,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,$migrtA,$migrtB,10,low,yes,$total10Lc\n");

		}
    }
    

    ### The below "until" statement (other half of the "do" statement earlier) keeps the simulation going until the ending conditions are met. 
    ### Here, it's when endgen has been reached 
    until($gen >= $endgen); 

print(FH "\n,propA1c,propA2c,propA3c,propA4c,propA5c,propA6c,propA7c,propA8c,propA9c,propA10c
            $simno,$propA1c,$propA2c,$propA3c,$propA4c,$propA5c,$propA6c,$propA7c,$propA8c,$propA9c,$propA10c\n");

print(FHd "\n,propA1c,propA2c,propA3c,propA4c,propA5c,propA6c,propA7c,propA8c,propA9c,propA10c
            $simno,$propA1c,$propA2c,$propA3c,$propA4c,$propA5c,$propA6c,$propA7c,$propA8c,$propA9c,$propA10c\n");

			#print(FH "simno,end1,end2,end3,end4,end5,end6,end7,end8,end9,end10,
            #$run,$end_1,$end_2,$end_3,$end_4,$end_5,$end_6,$end_7,$end_8,$end_9,$end_10\n"); 
 
   $finalgen = $gen; #the number of generations the simulation ran for 
    
   if($totalA == 0) { $extinctionA = 1; $total_extinctionsA++; }
   else { $extinctionA = 0; } 

   if($totalB == 0) { $extinctionB = 1; $total_extinctionsB++; }
   else { $extinctionB = 0; } 



    ### Print results to per-simulation file
    #print(FH1 "$total,$totalA,$totalB,$totalA1,$totalA2,$totalA3,$totalA4,$totalA5,$totalA6,$totalA7,$totalA8,$totalA9,$totalA10,$totalB1,$totalB2,$totalB3,$totalB4,$totalB5,$totalB6,$totalB7,$totalB8,$totalB9,$totalB10,$propATotal,$propBTotal,$finalgen,$simno\n");

    
    # accumulate results for summary file 
    $avg_total = $avg_total + $total; 
    $avg_total_A = $avg_total_A + $totalA; 
    $avg_total_B = $avg_total_B + $totalB; 
    $avg_prop_A = $avg_prop_A + $propATotal; 
    $avg_prop_B = $avg_prop_B + $propBTotal; 
    $avg_endgen = $avg_endgen + $finalgen; 
    
    $avgA1 = $avgA1 + $totalA1;
    $avgA2 = $avgA2 + $totalA2;
    $avgA3 = $avgA3 + $totalA3;
    $avgA4 = $avgA4 + $totalA4;
    $avgA5 = $avgA5 + $totalA5;
    $avgA6 = $avgA6 + $totalA6;
    $avgA7 = $avgA7 + $totalA7;
    $avgA8 = $avgA8 + $totalA8;
    $avgA9 = $avgA9 + $totalA9;
    $avgA10 = $avgA10 + $totalA10;
	
    
	
	
    # Find population sizes at the end generation
    if ($gen==$endgen){ # records population at the last generation
        $end_A1H = $total1Hv + $total1Hc; $end_A1M = $total1Mv + $total1Mc; $end_A1L = $total1Lv + $total1Lc;
        $end_A2H = $total2Hv + $total2Hc; $end_A2M = $total2Mv + $total2Mc; $end_A2L = $total2Lv + $total2Lc;
        $end_A3H = $total3Hv + $total3Hc; $end_A3M = $total3Mv + $total3Mc; $end_A3L = $total3Lv + $total3Lc;
        $end_A4H = $total4Hv + $total4Hc; $end_A4M = $total4Mv + $total4Mc; $end_A4L = $total4Lv + $total4Lc;
        $end_A5H = $total5Hv + $total5Hc; $end_A5M = $total5Mv + $total5Mc; $end_A5L = $total5Lv + $total5Lc;
        $end_A6H = $total6Hv + $total6Hc; $end_A6M = $total6Mv + $total6Mc; $end_A6L = $total6Lv + $total6Lc;
        $end_A7H = $total7Hv + $total7Hc; $end_A7M = $total7Mv + $total7Mc; $end_A7L = $total7Lv + $total7Lc;
        $end_A8H = $total8Hv + $total8Hc; $end_A8M = $total8Mv + $total8Mc; $end_A8L = $total8Lv + $total8Lc;
        $end_A9H = $total9Hv + $total9Hc; $end_A9M = $total9Mv + $total9Mc; $end_A9L = $total9Lv + $total9Lc;
        $end_A10H = $total10Hv + $total10Hc; $end_A10M = $total10Mv + $total10Mc; $end_A10L = $total10Lv + $total10Lc;

		
		# Total end populations, combining each dispersal phenotype
        $end_1 = $end_A1H + $end_A1M + $end_A1L;
        $end_2 = $end_A2H + $end_A2M + $end_A2L;
        $end_3 = $end_A3H + $end_A3M + $end_A3L;
        $end_4 = $end_A4H + $end_A4M + $end_A4L;
        $end_5 = $end_A5H + $end_A5M + $end_A5L;
        $end_6 = $end_A6H + $end_A6M + $end_A6L;
        $end_7 = $end_A7H + $end_A7M + $end_A7L;
        $end_8 = $end_A8H + $end_A8M + $end_A8L;
        $end_9 = $end_A9H + $end_A9M + $end_A9L;
        $end_10 = $end_A10H + $end_A10M + $end_A10L;
		
		
		# # Find allele frequencies for each simulation
		# if statements to avoid division by 0 in the case that end populations are 0

	if($end_1>100){
	$proc1=$propA1c+$proc1;}
	if($end_2>100){
	$proc2=$propA2c+$proc2;}
	if($end_3>100){
	$proc3=$propA3c+$proc3;}
	if($end_4>100){
	$proc4=$propA4c+$proc4;}
	if($end_5>100){
	$proc5=$propA5c+$proc5;}
	if($end_6>100){
	$proc6=$propA6c+$proc6;}
	if($end_7>100){
	$proc7=$propA7c+$proc7;}
	if($end_8>100){
	$proc8=$propA8c+$proc8;}
	if($end_9>100){
	$proc9=$propA9c+$proc9;}
	if($end_10>100){
	$proc10=$propA10c+$proc10;}
	
	if($end_1>100){
		$freq_1H = $end_A1H/$end_1; $freq_1M = $end_A1M/$end_1; $freq_1L = $end_A1L/$end_1;++$runs1; 
		}else{$freq_1H = $freq_1M = $freq_1L = 0;}	
        
		if($end_2>100){
		$freq_2H = $end_A2H/$end_2; $freq_2M = $end_A2M/$end_2; $freq_2L = $end_A2L/$end_2;++$runs2;
		}else{$freq_2H = $freq_2M = $freq_2L = 0;}
        
		if($end_3>100){
        $freq_3H = $end_A3H/$end_3; $freq_3M = $end_A3M/$end_3; $freq_3L = $end_A3L/$end_3;++$runs3;
		}else{$freq_3H = $freq_3M = $freq_3L = 0;}
		
		if($end_4>100){
        $freq_4H = $end_A4H/$end_4; $freq_4M = $end_A4M/$end_4; $freq_4L = $end_A4L/$end_4;++$runs4;
		}else{$freq_4H = $freq_4M = $freq_4L = 0;}
		
		if($end_5>100){
        $freq_5H = $end_A5H/$end_5; $freq_5M = $end_A5M/$end_5; $freq_5L = $end_A5L/$end_5;++$runs5;
		}else{$freq_5H = $freq_5M = $freq_5L = 0;}
		
		if($end_6>100){
        $freq_6H = $end_A6H/$end_6; $freq_6M = $end_A6M/$end_6; $freq_6L = $end_A6L/$end_6;++$runs6;
		}else{$freq_6H = $freq_6M = $freq_6L = 0;}
		
		if($end_7>100){
        $freq_7H = $end_A7H/$end_7; $freq_7M = $end_A7M/$end_7; $freq_7L = $end_A7L/$end_7;++$runs7;
		}else{$freq_7H = $freq_7M = $freq_7L = 0;}
		
		if($end_8>100){
        $freq_8H = $end_A8H/$end_8; $freq_8M = $end_A8M/$end_8; $freq_8L = $end_A8L/$end_8;++$runs8;
		}else{$freq_8H = $freq_8M = $freq_8L = 0;}
		
		if($end_9>100){
		 $freq_9H = $end_A9H/$end_9; $freq_9M = $end_A9M/$end_9; $freq_9L = $end_A9L/$end_9;++$runs9;
		}else{$freq_9H = $freq_9M = $freq_9L = 0;}
       
		if($end_10>100){
        $freq_10H = $end_A10H/$end_10; $freq_10M = $end_A10M/$end_10; $freq_10L = $end_A10L/$end_10;++$runs10;
		}else{$freq_10H = $freq_10M = $freq_10L = 0;}
	

        $avg_1H = $avg_1H + $freq_1H; $avg_1M = $avg_1M + $freq_1M; $avg_1L = $avg_1L + $freq_1L;
        $avg_2H = $avg_2H + $freq_2H; $avg_2M = $avg_2M + $freq_2M; $avg_2L = $avg_2L + $freq_2L;
        $avg_3H = $avg_3H + $freq_3H; $avg_3M = $avg_3M + $freq_3M; $avg_3L = $avg_3L + $freq_3L;
        $avg_4H = $avg_4H + $freq_4H; $avg_4M = $avg_4M + $freq_4M; $avg_4L = $avg_4L + $freq_4L;
        $avg_5H = $avg_5H + $freq_5H; $avg_5M = $avg_5M + $freq_5M; $avg_5L = $avg_5L + $freq_5L;
        $avg_6H = $avg_6H + $freq_6H; $avg_6M = $avg_6M + $freq_6M; $avg_6L = $avg_6L + $freq_6L;
        $avg_7H = $avg_7H + $freq_7H; $avg_7M = $avg_7M + $freq_7M; $avg_7L = $avg_7L + $freq_7L;
        $avg_8H = $avg_8H + $freq_8H; $avg_8M = $avg_8M + $freq_8M; $avg_8L = $avg_8L + $freq_8L;
        $avg_9H = $avg_9H + $freq_9H; $avg_9M = $avg_9M + $freq_9M; $avg_9L = $avg_9L + $freq_9L;
        $avg_10H = $avg_10H + $freq_10H; $avg_10M = $avg_10M + $freq_10M; $avg_10L = $avg_10L + $freq_10L;
print(FH "simno,end_A1H,end_A1M,end_A1L,end_A2H,end_A2M,end_A2L,end_A3H,end_A3M,end_A3L,end_A4H,end_A4M,end_A4L,end_A5H,end_A5M,end_A5L,end_A6H,end_A6M,end_A6L,end_A7H,end_A7M,end_A7L,end_A8H,end_A8M,end_A8L,end_A9H,end_A9M,end_A9L,end_A10H,end_A10M,end_A10L\n");

        print(FH "$simno,$end_A1H,$end_A1M,$end_A1L,$end_A2H,$end_A2M,$end_A2L,$end_A3H,$end_A3M,$end_A3L,$end_A4H,$end_A4M,$end_A4L,$end_A5H,$end_A5M,$end_A5L,$end_A6H,$end_A6M,$end_A6L,$end_A7H,$end_A7M,$end_A7L,$end_A8H,$end_A8M,$end_A8L,$end_A9H,$end_A9M,$end_A9L,$end_A10H,$end_A10M,$end_A10L");

#print(FH "\n,propA1c,propA2c,propA3c,propA4c,propA5c,propA6c,propA7c,propA8c,propA9c,propA10c
 #           $simno,$propA1c,$propA2c,$propA3c,$propA4c,$propA5c,$propA6c,$propA7c,$propA8c,$propA9c,$propA10c\n");
		print(FH "simno,end1,end2,end3,end4,end5,end6,end7,end8,end9,end10,\n");
 print(FH "simno,end1,end2,end3,end4,end5,end6,end7,end8,end9,end10,
            $simno,$end_1,$end_2,$end_3,$end_4,$end_5,$end_6,$end_7,$end_8,$end_9,$end_10\n");
print(FH5"$simno,$propA1c,$end_1\n,$propA2c,$end_2\n,$propA3c,$end_3\n,$propA4c,$end_4\n,$propA5c,$end_5\n,$propA6c,$end_6\n,$propA7c,$end_7\n,$propA8c,$end_8\n,$propA9c,$end_9\n,$propA10c,$end_10\n");

print(FH "simno,end_A1H,end_A1M,end_A1L,end_A2H,end_A2M,end_A2L,end_A3H,end_A3M,end_A3L,end_A4H,end_A4M,end_A4L,end_A5H,end_A5M,end_A5L,end_A6H,end_A6M,end_A6L,end_A7H,end_A7M,end_A7L,end_A8H,end_A8M,end_A8L,end_A9H,end_A9M,end_A9L,end_A10H,end_A10M,end_A10L\n");

        print(FHd "$simno,$end_A1H,$end_A1M,$end_A1L,$end_A2H,$end_A2M,$end_A2L,$end_A3H,$end_A3M,$end_A3L,$end_A4H,$end_A4M,$end_A4L,$end_A5H,$end_A5M,$end_A5L,$end_A6H,$end_A6M,$end_A6L,$end_A7H,$end_A7M,$end_A7L,$end_A8H,$end_A8M,$end_A8L,$end_A9H,$end_A9M,$end_A9L,$end_A10H,$end_A10M,$end_A10L");
		print(FHd "simno,end1,end2,end3,end4,end5,end6,end7,end8,end9,end10,\n");
 print(FHd "simno,end1,end2,end3,end4,end5,end6,end7,end8,end9,end10,
            $simno,$end_1,$end_2,$end_3,$end_4,$end_5,$end_6,$end_7,$end_8,$end_9,$end_10\n");
print(FH5d"$simno,$propA1c,$end_1\n,$propA2c,$end_2\n,$propA3c,$end_3\n,$propA4c,$end_4\n,$propA5c,$end_5\n,$propA6c,$end_6\n,$propA7c,$end_7\n,$propA8c,$end_8\n,$propA9c,$end_9\n,$propA10c,$end_10\n");


		}    


    ### Increase simulation number by 1 and go onto the next run of the simulation 
    $simno++;
    } # End of while-loop (the full simulation)


#calculate average statistics for summary file
$avg_total = $avg_total/$runs; 
$avg_total_A = $avg_total_A/$runs; 
$avg_total_B = $avg_total_B/$runs; 
$avg_prop_A = $avg_prop_A/$runs; 
$avg_prop_B = $avg_prop_B/$runs; 
$avg_endgen = $avg_endgen/$runs; 
		if($runs1>0){
$avgA1 = $avgA1/$runs1; $avgB1 = $avgB1/$runs;}
	if($runs2>0){
$avgA2 = $avgA2/$runs2; $avgB2 = $avgB2/$runs; }
	if($runs3>0){
$avgA3 = $avgA3/$runs3; $avgB3 = $avgB3/$runs; }
	if($runs4>0){
$avgA4 = $avgA4/$runs4; $avgB4 = $avgB4/$runs; }
	if($runs5>0){
$avgA5 = $avgA5/$runs5; $avgB5 = $avgB5/$runs; }
	if($runs6>0){
$avgA6 = $avgA6/$runs6; $avgB6 = $avgB6/$runs; }
	if($runs7>0){
$avgA7 = $avgA7/$runs7; $avgB7 = $avgB7/$runs; }
	if($runs8>0){
$avgA8 = $avgA8/$runs8; $avgB8 = $avgB8/$runs; }
	if($runs9>0){
$avgA9 = $avgA9/$runs9; $avgB9 = $avgB9/$runs; }
	if($runs10>0){
$avg10 = $avg10/$runs10; $avgB10 = $avgB10/$runs; }

	if($runs1>0){
$av_1H = $avg_1H/$runs1; $av_1M = $avg_1M/$runs1; $av_1L = $avg_1L/$runs1;++$ct;}
	if($runs2>0){
$av_2H = $avg_2H/$runs2; $av_2M = $avg_2M/$runs2; $av_2L = $avg_2L/$runs2;++$ct;}
	if($runs3>0){
$av_3H = $avg_3H/$runs3; $av_3M = $avg_3M/$runs3; $av_3L = $avg_3L/$runs3;++$ct;}
	if($runs4>0){
$av_4H = $avg_4H/$runs4; $av_4M = $avg_4M/$runs4; $av_4L = $avg_4L/$runs4;++$ct;}
	if($runs5>0){
$av_5H = $avg_5H/$runs5; $av_5M = $avg_5M/$runs5; $av_5L = $avg_5L/$runs5;++$ct;}
	if($runs6>0){
$av_6H = $avg_6H/$runs6; $av_6M = $avg_6M/$runs6; $av_6L = $avg_6L/$runs6;++$ct;}
	if($runs7>0){
$av_7H = $avg_7H/$runs7; $av_7M = $avg_7M/$runs7; $av_7L = $avg_7L/$runs7;++$ct;}
	if($runs8>0){
$av_8H = $avg_8H/$runs8; $av_8M = $avg_8M/$runs8; $av_8L = $avg_8L/$runs8;++$ct;}
	if($runs9>0){
$av_9H = $avg_9H/$runs9; $av_9M = $avg_9M/$runs9; $av_9L = $avg_9L/$runs9;++$ct;}
	if($runs10>0){
$av_10H = $avg_10H/$runs10; $av_10M = $avg_10M/$runs10; $av_10L = $avg_10L/$runs10;++$ct;}
$proc1av=$proc1/$runs1;
$proc2av=$proc2/$runs2;
$proc3av=$proc3/$runs3;
$proc4av=$proc4/$runs4;
$proc5av=$proc5/$runs5;
$proc6av=$proc6/$runs6;
$proc7av=$proc7/$runs7;
$proc8av=$proc8/$runs8;
$proc9av=$proc9/$runs9;
$proc10av=$proc10/$runs10;


$master_freq_H = ($av_1H + $av_2H + $av_3H + $av_4H + $av_5H + $av_6H + $av_7H + $av_8H + $av_9H + $av_10H) / $ct;
$master_freq_M = ($av_1M + $av_2M + $av_3M + $av_4M + $av_5M + $av_6M + $av_7M + $av_8M + $av_9M + $av_10M) / $ct;
$master_freq_L = ($av_1L + $av_2L + $av_3L + $av_4L + $av_5L + $av_6L + $av_7L + $av_8L + $av_9L + $av_10L) / $ct;

 
### Print cumulative results to summary file
print(FH "\nrun,avg_1H,avg_1M,avg_1L,avg_2H,avg_2M,avg_2L,avg_3H,avg_3M,avg_3L,avg_4H,avg_4M,avg_4L,avg_5H,avg_5M,avg_5L,avg_6H,avg_6M,avg_6L,avg_7H,avg_7M,avg_7L,avg_8H,avg_8M,avg_8L,avg_9H,avg_9M,avg_9L,avg_10H,avg_10M,avg_10L,freq_10H,freq_10M,freq_10L
            $runs,$av_1H,$av_1M,$av_1L,$av_2H,$av_2M,$av_2L,$av_3H,$av_3M,$av_3L,$av_4H,$av_4M,$av_4L,$av_5H,$av_5M,$av_5L,$av_6H,$av_6M,$av_6L,$av_7H,$av_7M,$av_7L,$av_8H,$av_8M,$av_8L,$av_9H,$av_9M,$av_9L,$av_10H,$av_10M,$av_10L\n\n
            avg. total pop size,total extinctions\n$avg_total_A,$total_extinctionsA\n");
#print "runs1=$runs1, runs2=$runs2,runs3=$runs3,runs4=$runs4";
#print "avg_1h=$avg_1H,avg_2h=$avg_2H;avg_3h=$avg_3H;avg_4h=$avg_4H";			
### Print cumulative results to summary file
print(FH "\nrun,proc1c,proc2c,proc3c,proc4c,proc5c,proc6c,proc7c,proc8c,proc9c,proc10c,
            $simno,$proc1av,$proc2av,$proc3av,$proc4av,$proc5av,$proc6av,$proc7av,$proc8av,$proc9av,$proc10av\n
            master_freq_H,master_freq_M,master_freq_L\n$master_freq_H,$master_freq_M,$master_freq_L");			
			
print(FHd "\nrun,avg_1H,avg_1M,avg_1L,avg_2H,avg_2M,avg_2L,avg_3H,avg_3M,avg_3L,avg_4H,avg_4M,avg_4L,avg_5H,avg_5M,avg_5L,avg_6H,avg_6M,avg_6L,avg_7H,avg_7M,avg_7L,avg_8H,avg_8M,avg_8L,avg_9H,avg_9M,avg_9L,avg_10H,avg_10M,avg_10L,freq_10H,freq_10M,freq_10L
            $runs,$av_1H,$av_1M,$av_1L,$av_2H,$av_2M,$av_2L,$av_3H,$av_3M,$av_3L,$av_4H,$av_4M,$av_4L,$av_5H,$av_5M,$av_5L,$av_6H,$av_6M,$av_6L,$av_7H,$av_7M,$av_7L,$av_8H,$av_8M,$av_8L,$av_9H,$av_9M,$av_9L,$av_10H,$av_10M,$av_10L\n\n
            avg. total pop size,total extinctions\n$avg_total_A,$total_extinctionsA\n");

print(FHd "\nrun,proc1c,proc2c,proc3c,proc4c,proc5c,proc6c,proc7c,proc8c,proc9c,proc10c,
            $simno,$proc1av,$proc2av,$proc3av,$proc4av,$proc5av,$proc6av,$proc7av,$proc8av,$proc9av,$proc10av\n
            master_freq_H,master_freq_M,master_freq_L\n$master_freq_H,$master_freq_M,$master_freq_L");			
	
			
			#print(FH "\n,propA1c,propA2c,propA3c,propA4c,propA5c,propA6c,propA7c,propA8c,propA9c,propA10c
 #           $simno,$propA1c,$propA2c,$propA3c,$propA4c,$propA5c,$propA6c,$propA7c,$propA8c,$propA9c,$propA10c\n");
	#		print(FH "simno,end1,end2,end3,end4,end5,end6,end7,end8,end9,end10,
     #       $simno,$end_1,$end_2,$end_3,$end_4,$end_5,$end_6,$end_7,$end_8,$end_9,$end_10\n");
#print(FH5"$simno,$propA1c,$end_1\n,$propA2c,$end_2\n,$propA3c,$end_3\n,$propA4c,$end_4\n,$propA5c,$end_5\n,$propA6c,$end_6\n,$propA7c,$end_7\n,$propA8c,$end_8\n,$propA9c,$end_9\n,$propA10c,$end_10\n");
	 
			
# End that particular set of parameters; move on to the next set of parameters (if looping, else simulation is over) 

}   # Closing bracket for file looping 
}

# } # Closing bracket for mutation rate loop
# } # For the migration rate loop
# } # For cancost 
# } # For num_can_offspring
# } # For carrying capacity 
# } # For dead 
