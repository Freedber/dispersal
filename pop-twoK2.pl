

# # # Two Species # # #

# Two species go head to head in this once in a lifetime competition to the death.  Who will win?  Species A
# is equiped with a higher migration rate, better to infultrate other habitats, but that migration rate could
# be its undoing.  Or will the more stationary Species B claim the day?  By the end of 
# this 20-minute simulation, who will be the victor?  Stay tuned to hear this story of love and loss, 
# bloodthirst and betrayal, infultration and intrigue whose ultimate resolution truly is a matter of life and 
# death.


############################################################
######################## Parameters ########################
############################################################
open(FH6, ">allruns.csv");
print (FH6 "cost,benefit,Total IC>0,Total IC>0.5,empty,empty patches A,empty patches B, Species A extinct, Species B extinct,original A absent, original B absent\n");
for ($test=0;$test<1;$test=$test+0.1){
for ($bene=0;$bene<1;$bene=$bene+0.1){
use constant RANDBITS => 15; 
use constant RAND_MAX => 2**RANDBITS;

srand(time ^ ($$ + ($$ << 15)) ); # Seeds the PRNG

### Logistics of model 
$runs = 100; # Number of runs
$invgen = 0; # Generation when individuals can start migrating
$endgen = 1000; 
#open(FH7, ">>G:/Shared drives/sim1119/test1.csv");
#print (FH7 "ttt");

### Model Parameters 
$carrying_capacity = 100; # Carrying capacity 
$mutrt = 0.000025; # Nuclear mutation rate   start with 10-6 but then switch to 10-5
#$mutrt = 0.000002; # Nuclear mutation rate   start with 10-6 but then switch to 10-5

#$cancost_rand = 6.05; #Used to seed the random number generator to determine the $cancost value each time the mother has a cannibal phenotype --cancost = 1 below 
#$dead = 0.25; #the probability of expressing the cannibal phenotype if an individual has the cannibal genotype
$num_cannibal_offspring = 1; #the number of additional offspring that a cannibal female can produce each mating 

$chance_cannibalism_cost = 0;  #probability that a cannibal mother will eat others of her population 
$chance_cannibalism_offspring = 0;  #probability that a cannibal mother will have additional offspring 


# Species A and B differ in migration rates 
$migrtA = 0.00;  # High Dispersal Migration Rate
$migrtB = 0.000;  # Low Dispersal Migration Rate
# open(FH2, ">posterT.csv"); 

### Reproductive Control 
# Used to determine the number of offspring each female will have (average offspring per female is 2.05 1 (2.2) )

$Rexp = 20; 
$master_R = (2.15); 

$femrep = 1; # Control of reproductive output; 1=females, 0=both sexes
$offSR = 0.5; # Average proportion of females in offspring 


### Fitness Costs 
# CURRENTLY NOT IN USE
$flatFC = 0.0; # Flat fitness cost (probability of offspring death when mother is a non-cannibal) 
$canFC0 = $canFC1 = $canFC2 = 0.0; # Cannibal fitness cost (probability of offspring death when cannibal mother is of non-cannibal genotype or cannibal genotype) 

$costs=$test;
$benefits=$bene;

@all_loop_costs = ($costs); # cost to the population from having a cannibalistic individual
@all_loop_benefits = ($benefits); # benefit to individual for being a cannibal

$extaltot=0;
$extalA=0;
$extalB=0;
$cantal50A=0;
$cantal50B=0;
$cantalA=0;
$cantalB=0;
$extorA=0;
$extorB=0;

# Immigration Rates--don't mess with 
$staySamePopA = 1 - $migrtA;
$moveNewPopA = (1 - $staySamePopA) / 9;

$staySamePopB = 1 - $migrtB;
$moveNewPopB = (1 - $staySamePopB) / 9;  



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

sub countCannibals(\@@){

    # This subroutine will return the number of cannibals in a given population array 
    
    # Arguments, in order: 
        # Population array for which a cannibal count will be generated 
        
    $pop = @_[0]; 
    
    $can_count = 0; 
    foreach(@$pop){
        if($$_[0][0][0] == 1 || $$_[0][0][1] == 1) {$can_count++; }
    }
    
    return $can_count; 
}

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
    if ($gender > $offSR) { push(@$pushMale, 0);} # Offspring is female
    else { push(@$pushFemale, 0); } # Offspring is male 
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
    
    ### Locus 0
    # Determine which chromosome (0/1) is inherited from each parent
    $ch0F = int(rand(2)); # From mother
    $ch0M = int(rand(2)); # From father
    # Send relevant allele to offspring nucleus array
    push(@{$offN[0][0]}, $$mother[0][0][$ch0F]);
    push(@{$offN[0][0]}, $$father[0][0][$ch0M]);
    
    
    ##################################################
    ########## Potential nuclear mutation ############
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
        # species code (e.g, 'A')
    
    $origin = $_[0];
    $destination = $_[1];
    $species_code = $_[2]; 
    
    if($species_code eq 'A') {
        if($origin == $destination) {$immRate = $staySamePopA;} 
        else {$immRate = $moveNewPopA;}
    }
    else{
        if($origin == $destination) {$immRate = $staySamePopB;}
        else {$immRate = $moveNewPopB; }
    }
    
    return $immRate;
}


sub migration ($$\@\@\@\@\@\@\@\@\@\@\@\@\@\@\@\@\@\@\@\@\@\@\@\@\@) {
    # Used to randomly move individuals from one habitat to another.
    
    # Arguments, in order:
        # Habitat number (e.g., 1)
        # Species code (e.g., 'A') 
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
    $species_code = $_[1]; 
    
    $mothers = @_[2]; 
    $female_tracker = @_[3]; 
    $male_tracker = @_[4]; 
       
    $push_FemH1 = @_[5]; 
    $push_FemH2 = @_[6]; 
    $push_FemH3 = @_[7]; 
    $push_FemH4 = @_[8]; 
    $push_FemH5 = @_[9]; 
    $push_FemH6 = @_[10]; 
    $push_FemH7 = @_[11]; 
    $push_FemH8 = @_[12]; 
    $push_FemH9 = @_[13]; 
    $push_FemH10 = @_[14]; 
    
    $push_MalH1 = @_[15]; 
    $push_MalH2 = @_[16]; 
    $push_MalH3 = @_[17]; 
    $push_MalH4 = @_[18]; 
    $push_MalH5 = @_[19]; 
    $push_MalH6 = @_[20]; 
    $push_MalH7 = @_[21]; 
    $push_MalH8 = @_[22]; 
    $push_MalH9 = @_[23]; 
    $push_MalH10 = @_[24]; 
    
    $femaleOffspring = @_[25]; 
    $maleOffspring = @_[26]; 
    
    $num_mothers = @$mothers; 
    
    
    ### Based on immigration rates, set probabilites for sending an individual to each habitat
    $imm1 = immigrationRates($habitatNumber, 1, $species_code);
    $imm2 = $imm1 + immigrationRates($habitatNumber, 2, $species_code);
    $imm3 = $imm2 + immigrationRates($habitatNumber, 3, $species_code);
    $imm4 = $imm3 + immigrationRates($habitatNumber, 4, $species_code);
    $imm5 = $imm4 + immigrationRates($habitatNumber, 5, $species_code);
    $imm6 = $imm5 + immigrationRates($habitatNumber, 6, $species_code);
    $imm7 = $imm6 + immigrationRates($habitatNumber, 7, $species_code);
    $imm8 = $imm7 + immigrationRates($habitatNumber, 8, $species_code); 
    $imm9 = $imm8 + immigrationRates($habitatNumber, 9, $species_code); 
    
    
    $current_female = 0; $current_male = 0; 

    ### Migrate each individual to new array (which could still be native habitat) 
    for($mother = 0; $mother < $num_mothers; $mother++) { 
        $migrate = double_rand();
        
        # check the tracker numbers for this mother        
        $num_female_children = $female_tracker->[$mother];  
        $num_male_children = $male_tracker->[$mother]; 
        #print "children: $num_female_children and $num_male_children\n"; 

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

sub makeNextGenVegan($$$\@\@){
    #Creates the next generation of the Vegan species for a given habitat's population of Vegans.  
    #Irrelevant in this scenario because there are no Vegans. 
    
    # Arguments, in order: 
        # Number of females (e.g., $lenFV1)
        # Number of males (e.g., $lenMV1)
        # Maximum offspring (e.g., $max_offspringV1)
        # Female offspring array (e.g., @FoffspringV1)
        # Male offspring array (e.g., @MoffspringV1)
        
        $lenFem = @_[0]; 
        $lenMal = @_[1]; 
        $max_off = @_[2];
        $female_offspring = @_[3]; 
        $male_offspring = @_[4]; 
 
    
    ### Vegan species in Habitat X is invariably non-cannibal; thus, don't need to randomly choose parents  
        if($lenFem > 0 && $lenMal > 0) {
            for($offs = 0; $offs < $max_off; $offs++) { # -0.5 is for rounding correction
                createVeganOffspring(@$female_offspring, @$male_offspring); 
            }
        }    
}

sub makeNextGenCarnivore($$$\@\@\@\@) {
    #Creates the next generation of the cannibal species for a given habitat's population of cannibals. 
    
    # Arguments, in order: 
        # Number of females (e.g., $lenFC1)
        # Number of males (e.g., $lenMC1)
        # Maximum offspring (e.g., $max_offspringC1)
        # Female adults array (e.g., @femalesC1)
        # Male adults array (e.g., @malesC1)
        # Female offspring array (e.g., @FoffspringC1)
        # Male offspring array (e.g., @MoffspringC1)
           
        $lenFem = $_[0]; 
        $lenMal = $_[1]; 
        $max_off = $_[2]; 
        $females = @_[3]; 
        $males = @_[4]; 
        $female_offspring = @_[5]; 
        $male_offspring = @_[6]; 
        
        $offs = 0; #number of offspring produced during mating season    
        if($lenFem > 0 && $lenMal > 0) {
            for($fem = 0; $fem < $lenFem; $fem++) {
                
                if($offs >= $carrying_capacity) { $fem = $lenFem; $offs2 = $R; $add_off = $num_cannibal_offspring; $die = $cancost; } 
                #create first offspring
                if($lenFem > 0 && $lenMal > 0){
                    $Fmate = int(rand($lenFem));
                    $Mmate = int(rand($lenMal));
            }
                
                #the reproductive rate for this particular female; represents the number of offspring she will have    
 
                $size = $lenFem + $lenMal; 
                $R = calcR($size);              
           
                 for($offs2 = 0; $offs2 < $R; $offs2++){
                    if($offs >= $carrying_capacity) { $fem = $lenFem; $offs2 = $R; $add_off = $num_cannibal_offspring; $die = $cancost; } 

                    createCarnivoreOffspring(@{$$females[$Fmate]}, @{$$males[$Mmate]}, @$female_offspring, @$male_offspring);
                    $offs = $offs + 1; 
                 }
            
                ### Does the mother in question possess a cannibal genotype? 1=yes; 0=no
                if($$females[$Fmate][0][0][0] == 1 || $$females[$Fmate][0][0][1] == 1) {$maternalCan = 1;} else {$maternalCan = 0;}

                # If the mother expresses the cannibal phenotype... 
                $ChanceCanOff = rand(); # random number between 0-1
                $ChanceCanCost = rand(); 
                $cancost = 1;   #$cancost = rand($cancost_rand);

                if($offs >= $carrying_capacity) { $fem = $lenFem; $offs2 = $R; $add_off = $num_cannibal_offspring; $die = $cancost; } 


                if($maternalCan == 1 && $ChanceCanOff < $chance_cannibalism_offspring){ #if cannibal mother can produce additional offspring  
                    # Cannibal mother's additional offspring
                    for($add_off = 0; $add_off < $num_cannibal_offspring; $add_off++){ #make the additional offsprnig that a cannibal produces
                        if($offs >= $carrying_capacity) { $fem = $lenFem; $offs2 = $R; $add_off = $num_cannibal_offspring; $die = $cancost; } 
 
                        createCarnivoreOffspring(@{$$females[$Fmate]}, @{$$males[$Mmate]}, @$female_offspring, @$male_offspring); #Second Offspring if the above is true 
                        $offs = $offs + 1; #counting additional offspring 
                    }
                } 
                        
                if($maternalCan == 1 && $ChanceCanCost < $chance_cannibalism_cost){ # if cannibal mother kills off others in her population 
                    # Individuals killed by cannibal mother when gathering resources to lay eggs 
                    for($die = 0; $die < $cancost; $die++){  
                        if($offs >= $carrying_capacity) { $fem = $lenFem; $offs2 = $R; $add_off = $num_cannibal_offspring; $die = $cancost; } 
   
                        $gender = int(rand(2));
                        if($gender == 0){$lenFem = $lenFem - 1; #she can kill off either males or females, and they are no longer able to mate (because they are dead) 
                            }else{$lenMal = $lenMal - 1;}
                    }                                
                }#end if the mother expresses the cannibal phenotype  
            

                if($lenMal <= 0 || $lenFem <= 0){$offs = $max_off; } #if the cannibals have killed off all the adults of one sex, then no more offspring can be produced                   
        }       
    }
}            

sub makeNextGenCompetingSpecies($$$$$\@\@\@\@\@\@\@\@\@\@\@\@){
    
    # This subroutine goes through the procedure for a given generation's mating season. 
    # It takes care of both Species A and Species B in the same habitat with a shared carrying capacity. 
    
    # Arguments, in order: 
        # Number of Species A females (e.g., $lenFA1)
        # Number of Species A males (e.g., $lenMA1)
        # Number of Species B females (e.g., $lenFB1)
        # Number of Species B males (e.g., $lenMB1)
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
           
        $lenFemA = $_[0];
        $lenMalA = $_[1];  
        $lenFemB = $_[2]; 
        $lenMalB = $_[3]; 
        $max_off = $_[4]; 
        $femalesA = @_[5]; 
        $malesA = @_[6];
        $femalesB = @_[7]; 
        $malesB = @_[8];  
        $female_offspringA = @_[9]; 
        $male_offspringA = @_[10]; 
        $female_offspringB = @_[11]; 
        $male_offspringB = @_[12]; 
        $female_trackerA = @_[13]; 
        $male_trackerA = @_[14]; 
        $female_trackerB = @_[15]; 
        $male_trackerB = @_[16]; 
        
        
        
        $females_remaining_A = $lenFemA;  #number of females left to produce offspring 
        $females_remaining_B = $lenFemB; 
        $off = 0; 
        
        # mating will continue until either the carrying capcity is reached or until all females from both species have produced offspring 
        while($off < $max_off && ($females_remaining_A > 0 || $females_remaining_B > 0)){
            
            #randomly choose which species' turn it is to reproduce 
            $which_species = rand(); 
            $sizeA = $lenFemA + $lenMalA; 
            $sizeB = $lenFemB + $lenMalB; 
            $total = $sizeA + $sizeB; 
            if($which_species < $sizeA/$total){ ### Species A reproduces 
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
            else{ ### Species B reproduces 
                if($lenFemB == 0 || $lenMalB == 0) {$females_remaining_B = 0; }
                if($females_remaining_B > 0){ #Species B can reproduce unless all their females have already produced offspring (or were eaten) 
                    
                    # produce a Species B offspring 
                    
                    #find the female and male mates 
                    $Fmate = int(rand($lenFemB));
                    $Mmate = int(rand($lenMalB));
                    
                    #calculate reproductive rate of that female
                    $size = $lenFemB + $lenMalB;  
                    $total = $size + $lenFemA + $lenMalA; 
                    $R = calcR($size, $total); 
                    
                    $male_countB = 0; $female_countB = 0; 
                    #make offspring (based on that reproductive rate) whether or not Species B female is cannibal or non-cannibal 
                    for($off3 = 0; $off3 < $R; $off3++){
                        $track_gender = createCarnivoreOffspring(@{$$femalesB[$Fmate]}, @{$$malesB[$Mmate]}, @$female_offspringB, @$male_offspringB); 
                        $off++; 
                        if($track_gender == 1) {$male_countB++; }
                        else {$female_countB++; }
                    }
                    
                    ### Does the mother in question possess a cannibal genotype? 1=yes; 0=no
                    if($$femalesB[$Fmate][0][0][0] == 1 || $$femalesB[$Fmate][0][0][1] == 1) {$maternalCan = 1;} else {$maternalCan = 0;}
                    
                    #if female is cannibal, then she will make additional offspring 
                    $ChanceCanOff = rand(); 
                    $ChanceCanCost = rand(); 
                    $cancost = 1; #$cancost = rand($cancost_rand);
                    
                    # If the mother expresses the cannibal phenotype... 
                    if($maternalCan == 1 && $ChanceCanOff < $chance_cannibalism_offspring){ #if the cannibal mother produces additional offspring 
                        # Cannibal mother's additional offspring
                        for($add_off = 0; $add_off < $num_cannibal_offspring; $add_off++){ #make the additional offsprnig that a cannibal produces
                            $track_gender = createCarnivoreOffspring(@{$$femalesB[$Fmate]}, @{$$malesB[$Mmate]}, @$female_offspringB, @$male_offspringB); #Second Offspring if the above is true 
                            $off++;
                            if($track_gender == 1) {$male_countB++; }
                            else {$female_countB++; }  
                        }
                    }                    
                    
                    if($maternalCan == 1 && $ChanceCanCost < $chance_cannibalism_cost){ #if the cannibal mother kills off others of her population 
                        # Individuals killed by cannibal mother when gathering resources to lay eggs 
                        for($die = 0; $die < $cancost; $die++){  
                            $gender = int(rand(2));
                            if($gender == 0){$lenFemB = $lenFemB - 1; $females_remaining_B -= 1; #she can kill off either males or females, and they are no longer able to mate (because they are dead) 
                                }else{$lenMalB = $lenMalB - 1;}
                        }                
                    }  
                        
                    if($lenMalB <= 0 || $lenFemB <= 0){$females_remaining_B = 0; }

                    $females_remaining_B -= 1;   
                    push(@$female_trackerB, $female_countB); 
                    push(@$male_trackerB, $male_countB);  
                }             
            }
        }
}    
$filename = "$carrying_capacity$migrtA$migrtB$costs$benefits";
$sumname = "sum_$filename";

############################################################
####################### File Setup #########################
############################################################

### Create files to save output

    # Simulation summary (multiple trials)- used with parameter loops
    open(FHd, ">$sumname.$endgen.summary.csv");
    # Per-simulation results
    open(FH1d, ">$sumname.$endgen.eachSim.csv");
	select((select(FH1d), $|=1)[0]);
    # Per generation results
    #open(FH2, ">scn5v3_gens.csv"); 
		open(FH5d, ">$sumname.$endgen.columns.csv");
		open(FH, ">G:/Shared drives/sim1119/twospecies/$sumname.$endgen.summary.csv");
		open(FH1, ">G:/Shared drives/sim1119/twospecies/$sumname.$endgen.eachSim.csv");
		open(FH5, ">G:/Shared drives/sim1119/twospecies/$sumname.$endgen.columns.csv");
	select((select(FH1), $|=1)[0]);		
		select((select(FH5), $|=1)[0]);
		print(FH5"simno,prop can A,total pop A,prop can B,total pop B\n");
	print(FH5d"simno,prop can A,total pop A,prop can B,total pop B\n");


### Print non-looping parameters and column descriptions to files
# Simulation summary results (use skip __ lines when reading in file, e.g. to R)
print(FH "### PARAMETERS ###\nScenario number = 5 (interspecific competition/metapopulation/migration)\nBeginning of migration = generation $invgen\nEnding generation = $endgen\nReproductive output = $femrep (control by: 1=females / 0=both sexes)\nOffspring sex ratio = $offSR\nFC scaling = $flatFC/$canFC0/$canFC1/$canFC2 (scalars)\n\ncarrying capacity = $carrying_capacity\nmutation rate = $mutrt\nmigration rate A = $migrtA\nmigration rate B = $migrtB\ncancost = $cancost_rand\ncan_off = $num_cannibal_offspring\ndead = $dead\n\n");
print(FH "### COLUMNS ###\navg_total = Average final individuals\navg_totalA = Average final Species A individuals\navg_totalB = Average final Species B individuals\navgA/B# = Average number of Species A or B in Habitat #\navg_propA = average proportion of Species A in metacommunity\navg_propB = average proportion of Species B in metacommunity\navg_endgen = the average ending generation among runs\nruns = Number of simulations run with given parameter set\n\n### END HEADER ###\n\n");
print(FHd "### PARAMETERS ###\nScenario number = 5 (interspecific competition/metapopulation/migration)\nBeginning of migration = generation $invgen\nEnding generation = $endgen\nReproductive output = $femrep (control by: 1=females / 0=both sexes)\nOffspring sex ratio = $offSR\nFC scaling = $flatFC/$canFC0/$canFC1/$canFC2 (scalars)\n\ncarrying capacity = $carrying_capacity\nmutation rate = $mutrt\nmigration rate A = $migrtA\nmigration rate B = $migrtB\ncancost = $cancost_rand\ncan_off = $num_cannibal_offspring\ndead = $dead\n\n");
print(FHd "### COLUMNS ###\navg_total = Average final individuals\navg_totalA = Average final Species A individuals\navg_totalB = Average final Species B individuals\navgA/B# = Average number of Species A or B in Habitat #\navg_propA = average proportion of Species A in metacommunity\navg_propB = average proportion of Species B in metacommunity\navg_endgen = the average ending generation among runs\nruns = Number of simulations run with given parameter set\n\n### END HEADER ###\n\n");

# Per-simulation results (skip=_)
print(FH1 "### PARAMETERS ###\nScenario number = 5 (interspecific competition/metapopulation/migration)\nBeginning of migration = generation $invgen\nEnding generation = $endgen\nReproductive output = $femrep (control by: 1=females / 0=both sexes)\nOffspring sex ratio = $offSR\nFC scaling = $flatFC/$canFC0/$canFC1/$canFC2 (scalars)\n\ncarrying capacity = $carrying_capacity\nmutation rate = $mutrt\nmigration rate A = $migrtA\nmigration rate B = $migrtB\ncancost = $cancost_rand\ncan_off = $num_cannibal_offspring\ndead = $dead\n\n");
print(FH1 "### COLUMNS ###\ntotal = # individuals when simulation ends\ntotalA = # of Species A when simulation ends\ntotalB = # of Species B when simulation ends\ntotalA# = number of Species A in Habitat #\ntotalB# = number of Species B in Habitat #\npropA = proportion of Species A in metacommunity\npropB = proportion of Species B in metacommunity\nendgen = Final generation\nrun = Run number with that set of parameters\n\n### END HEADER ###\n\n");
print(FH1d "### PARAMETERS ###\nScenario number = 5 (interspecific competition/metapopulation/migration)\nBeginning of migration = generation $invgen\nEnding generation = $endgen\nReproductive output = $femrep (control by: 1=females / 0=both sexes)\nOffspring sex ratio = $offSR\nFC scaling = $flatFC/$canFC0/$canFC1/$canFC2 (scalars)\n\ncarrying capacity = $carrying_capacity\nmutation rate = $mutrt\nmigration rate A = $migrtA\nmigration rate B = $migrtB\ncancost = $cancost_rand\ncan_off = $num_cannibal_offspring\ndead = $dead\n\n");
print(FH1d "### COLUMNS ###\ntotal = # individuals when simulation ends\ntotalA = # of Species A when simulation ends\ntotalB = # of Species B when simulation ends\ntotalA# = number of Species A in Habitat #\ntotalB# = number of Species B in Habitat #\npropA = proportion of Species A in metacommunity\npropB = proportion of Species B in metacommunity\nendgen = Final generation\nrun = Run number with that set of parameters\n\n### END HEADER ###\n\n");

# Per generation results (skip=_)
print(FH2 "### PARAMETERS ###\nScenario number = 5 (interspecific competition/metapopulation/migration)\nBeginning of migration = generation $invgen\nEnding generation = $endgen\nReproductive output = $femrep (control by: 1=females / 0=both sexes)\nOffspring sex ratio = $offSR\nFC scaling = $flatFC/$canFC0/$canFC1/$canFC2 (scalars)\n\ncarrying capacity = $carrying_capacity\nmutation rate = $mutrt\nmigration rate A = $migrtA\nmigration rate B = $migrtB\ncancost = $cancost_rand\ncan_off = $num_cannibal_offspring\ndead = $dead\n\n");
print(FH2 "### COLUMNS ###\nTotalA/B#c/n = number of cannibal/non-cannibal individuals of Species A or B in Habitat\navgA/B# = avg number of individuals of Species A or B in Habitat #\ngen = generation #\nrun = Run number of the simulation\n\n### END HEADER ###\n\n");

# Per generation results (skip=_)
print(FH3 "### PARAMETERS ###\nScenario number = 5 (interspecific competition/metapopulation/migration)\nBeginning of migration = generation $invgen\nEnding generation = $endgen\nReproductive output = $femrep (control by: 1=females / 0=both sexes)\nOffspring sex ratio = $offSR\nFC scaling = $flatFC/$canFC0/$canFC1/$canFC2 (scalars)\n\ncarrying capacity = $carrying_capacity\nmutation rate = $mutrt\nmigration rate A = $migrtA\nmigration rate B = $migrtB\ncancost = $cancost_rand\ncan_off = $num_cannibal_offspring\ndead = $dead\n\n");
print(FH3 "### COLUMNS ###\nTotalA/B#c/n = number of cannibal/non-cannibal individuals of Species A or B in Habitat\navgA/B# = avg number of individuals of Species A or B in Habitat #\ngen = generation #\nrun = Run number of the simulation\n\n### END HEADER ###\n\n");


### Print column name headers to files
print(FH "avg_total,avg_totalA,avg_totalB,avgA1,avgA2,avgA3,avgA4,avgA5,avgA6,avgA7,avgA8,avgA9,avgA10,avgB1,avgB2,avgB3,avgB4,avgB5,avgB6,avgB7,avgB8,avgB9,avgB10,avg_propA,avg_propB,avg_endgen,runs\n");
print(FH1 "total,totalA,totalB,totalA1,totalA2,totalA3,totalA4,totalA5,totalA6,totalA7,totalA8,totalA9,totalA10,totalB1,totalB2,totalB3,totalB4,totalB5,totalB6,totalB7,totalB8,totalB9,totalB10,propA,propB,endgen,run\n");
print(FHd "avg_total,avg_totalA,avg_totalB,avgA1,avgA2,avgA3,avgA4,avgA5,avgA6,avgA7,avgA8,avgA9,avgA10,avgB1,avgB2,avgB3,avgB4,avgB5,avgB6,avgB7,avgB8,avgB9,avgB10,avg_propA,avg_propB,avg_endgen,runs,Acavg,Bcavg,ICpresent,Exttotal\n");
print(FH1d "total,totalA,totalB,totalA1,totalA2,totalA3,totalA4,totalA5,totalA6,totalA7,totalA8,totalA9,totalA10,totalB1,totalB2,totalB3,totalB4,totalB5,totalB6,totalB7,totalB8,totalB9,totalB10,ICA1,ICA2,ICA3,ICA4,ICA5,ICA6,ICA7,ICA8,ICA9,ICA10,ICB1,ICB2,ICB3,ICB4,ICB5,ICB6,ICB7,ICB8,ICB9,ICB10,propA,propB,endgen,run\n");

print(FH2 "totalA1c,totalA1n,totalA2c,totalA2n,totalA3c,totalA3n,totalA4c,totalA4n,totalA5c,totalA5n,totalA6c,totalA6n,totalA7c,totalA7n,totalA8c,totalA8n,totalA9c,totalA9n,totalA10c,totalA10n,totalB1c,totalB1n,totalB2c,totalB2n,totalB3c,totalB3n,totalB4c,totalB4n,totalB5c,totalB5n,totalB6c,totalB6n,totalB7c,totalB7n,totalB8c,totalB8n,totalB9c,totalB9n,totalB10c,totalB10n,gen,run\n"); 
print(FH3 "totalA1c,totalA1n,totalA2c,totalA2n,totalA3c,totalA3n,totalA4c,totalA4n,totalA5c,totalA5n,totalA6c,totalA6n,totalA7c,totalA7n,totalA8c,totalA8n,totalA9c,totalA9n,totalA10c,totalA10n,totalB1c,totalB1n,totalB2c,totalB2n,totalB3c,totalB3n,totalB4c,totalB4n,totalB5c,totalB5n,totalB6c,totalB6n,totalB7c,totalB7n,totalB8c,totalB8n,totalB9c,totalB9n,totalB10c,totalB10n,gen,run\n"); 





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

open(FH4d, ">$migrtA$migrtB.$endgen.csv");
open(FH4, ">G:/Shared drives/sim1119/twospecies/$migrtA$migrtB.$endgen.csv");
		
select((select(FH4), $|=1)[0]);
print(FH4 "gen,run,cost,benefit,species,habitat,cannibal,pop,\n");    
print(FH4d "gen,run,cost,benefit,species,habitat,cannibal,pop,\n");        
    
### Population survival indicators
$extinctionA = $extinctionB = 0; # simulation ends due to extinction of either Species A or Species B 
$total_extinctionsA = $total_extinctionsB = 0; # number of times Species A or Species B went extinct 

$simno = 1;

$avg_total = $avg_total_A = $avg_total_B = $avg_prop_A = $avg_prop_B = $avg_endgen = 0; 
$avgA1 = $avgA2 = $avgA3 = $avgA4 = $avgA5 = $avgA6 = $avgA7 = $avgA8 = $avgA9 = $avgB1 = $avgB2 = $avgB3 = $avgB4 = $avgB5 = $avgB6 = $avgB7 = $avgB8 = $avgB9 = 0; 


### Run through desired number of simulations
while ($simno <= $runs) {
    #print "File: $file, Simulation $simno:\n"; 
    print "Simulation $simno:\tMutation $mutrt; Migration A $migrtA; Migration B $migrtB; Cost $chance_cannibalism_cost; Offspring $chance_cannibalism_offspring; R exponent $Rexp; Capacity $carrying_capacity\n"; 
  


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
    
    # Male and female Species B populations of Habitats 1-9
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
    ### Species A non-cannibals: Habitas 1-5
    ### Species B non-cannibals: Habitats 6-10
    
    
    # Habitat 1
    for($ind = 0; $ind < $out; $ind++){
        if($mutrt > 0){
            $mutate0 = int(double_rand((1/$mutrt))); 
            $mutate1 = int(double_rand((1/$mutrt)));
        } 
        else { $mutate0 = 0; $mutate1 = 0; }
        if($mutate0 == 42 && $mutate1 == 42){
            push(@femalesA1, [ [[1, 1], [0, 0]], [[0]] ]); }
        elsif($mutate0 == 42 || $mutate1 == 42){
            push(@femalesA1, [ [[1, 0], [0, 0]], [[0]] ]); }            
        else{
            push(@femalesA1, [ [[0, 0], [0, 0]], [[0]] ]); } 
    }
    for($ind = 0; $ind < $k-$out; $ind++){ 
        if($mutrt > 0){
            $mutate0 = int(double_rand((1/$mutrt))); 
            $mutate1 = int(double_rand((1/$mutrt)));
        } 
        else { $mutate0 = 0; $mutate1 = 0; }
        if($mutate0 == 42 && $mutate1 == 42){
            push(@malesA1, [ [[1, 1], [0, 0]], [[0]] ]); }
        elsif($mutate0 == 42 || $mutate1 == 42){
            push(@malesA1, [ [[1, 0], [0, 0]], [[0]] ]); }        
        else{
            push(@malesA1, [ [[0, 0], [0, 0]], [[0]] ]); } 
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
            push(@femalesA3, [ [[1, 1], [0, 0]], [[0]] ]); }
        elsif($mutate0 == 42 || $mutate1 == 42){
            push(@femalesA3, [ [[1, 0], [0, 0]], [[0]] ]); }            
        else{
            push(@femalesA3, [ [[0, 0], [0, 0]], [[0]] ]); } 
    }
    for($ind = 0; $ind < $k-$out; $ind++){ 
        if($mutrt > 0){
            $mutate0 = int(double_rand((1/$mutrt))); 
            $mutate1 = int(double_rand((1/$mutrt)));
        } 
        else { $mutate0 = 0; $mutate1 = 0; }
        if($mutate0 == 42 && $mutate1 == 42){
            push(@malesA3, [ [[1, 1], [0, 0]], [[0]] ]); }
        elsif($mutate0 == 42 || $mutate1 == 42){
            push(@malesA3, [ [[1, 0], [0, 0]], [[0]] ]); }        
        else{
            push(@malesA3, [ [[0, 0], [0, 0]], [[0]] ]); } 
    }

    ### Habitat 4
    for($ind = 0; $ind < $out; $ind++){ 
        if($mutrt > 0){
            $mutate0 = int(double_rand((1/$mutrt))); 
            $mutate1 = int(double_rand((1/$mutrt)));
        } 
        else { $mutate0 = 0; $mutate1 = 0; }
        if($mutate0 == 42 && $mutate1 == 42){
            push(@femalesA4, [ [[1, 1], [0, 0]], [[0]] ]); }
        elsif($mutate0 == 42 || $mutate1 == 42){
            push(@femalesA4, [ [[1, 0], [0, 0]], [[0]] ]); }            
        else{
            push(@femalesA4, [ [[0, 0], [0, 0]], [[0]] ]); } 
    }
    for($ind = 0; $ind < $k-$out; $ind++){ 
        if($mutrt > 0){
            $mutate0 = int(double_rand((1/$mutrt))); 
            $mutate1 = int(double_rand((1/$mutrt)));
        } 
        else { $mutate0 = 0; $mutate1 = 0; }
        if($mutate0 == 42 && $mutate1 == 42){
            push(@malesA4, [ [[1, 1], [0, 0]], [[0]] ]); }
        elsif($mutate0 == 42 || $mutate1 == 42){
            push(@malesA4, [ [[1, 0], [0, 0]], [[0]] ]); }        
        else{
            push(@malesA4, [ [[0, 0], [0, 0]], [[0]] ]); } 
    }

    ### Habitat 5
    for($ind = 0; $ind < $out; $ind++){ 
        if($mutrt > 0){
            $mutate0 = int(double_rand((1/$mutrt))); 
            $mutate1 = int(double_rand((1/$mutrt)));
        } 
        else { $mutate0 = 0; $mutate1 = 0; } 
        if($mutate0 == 42 && $mutate1 == 42){
            push(@femalesA5, [ [[1, 1], [0, 0]], [[0]] ]); }
        elsif($mutate0 == 42 || $mutate1 == 42){
            push(@femalesA5, [ [[1, 0], [0, 0]], [[0]] ]); }            
        else{
            push(@femalesA5, [ [[0, 0], [0, 0]], [[0]] ]); } 
    }
    for($ind = 0; $ind < $k-$out; $ind++){ 
        if($mutrt > 0){
            $mutate0 = int(double_rand((1/$mutrt))); 
            $mutate1 = int(double_rand((1/$mutrt)));
        } 
        else { $mutate0 = 0; $mutate1 = 0; }
        if($mutate0 == 42 && $mutate1 == 42){
            push(@malesA5, [ [[1, 1], [0, 0]], [[0]] ]); }
        elsif($mutate0 == 42 || $mutate1 == 42){
            push(@malesA5, [ [[1, 0], [0, 0]], [[0]] ]); }        
        else{
            push(@malesA5, [ [[0, 0], [0, 0]], [[0]] ]); } 
    }

    ## Habitat 6
    for($ind = 0; $ind < $out; $ind++){ 
        if($mutrt > 0){
            $mutate0 = int(double_rand((1/$mutrt))); 
            $mutate1 = int(double_rand((1/$mutrt)));
        } 
        else { $mutate0 = 0; $mutate1 = 0; } 
        if($mutate0 == 42 && $mutate1 == 42){
            push(@femalesB6, [ [[1, 1], [0, 0]], [[0]] ]); }
        elsif($mutate0 == 42 || $mutate1 == 42){
            push(@femalesB6, [ [[1, 0], [0, 0]], [[0]] ]); }            
        else{
            push(@femalesB6, [ [[0, 0], [0, 0]], [[0]] ]); } 
    }
    for($ind = 0; $ind < $k-$out; $ind++){ 
        if($mutrt > 0){
            $mutate0 = int(double_rand((1/$mutrt))); 
            $mutate1 = int(double_rand((1/$mutrt)));
        } 
        else { $mutate0 = 0; $mutate1 = 0; }
        if($mutate0 == 42 && $mutate1 == 42){
            push(@malesB6, [ [[1, 1], [0, 0]], [[0]] ]); }
        elsif($mutate0 == 42 || $mutate1 == 42){
            push(@malesB6, [ [[1, 0], [0, 0]], [[0]] ]); }        
        else{
            push(@malesB6, [ [[0, 0], [0, 0]], [[0]] ]); } 
    }

    ## Habitat 7
    for($ind = 0; $ind < $out; $ind++){ 
        if($mutrt > 0){
            $mutate0 = int(double_rand((1/$mutrt))); 
            $mutate1 = int(double_rand((1/$mutrt)));
        } 
        else { $mutate0 = 0; $mutate1 = 0; }
        if($mutate0 == 42 && $mutate1 == 42){
            push(@femalesB7, [ [[1, 1], [0, 0]], [[0]] ]); }
        elsif($mutate0 == 42 || $mutate1 == 42){
            push(@femalesB7, [ [[1, 0], [0, 0]], [[0]] ]); }            
        else{
            push(@femalesB7, [ [[0, 0], [0, 0]], [[0]] ]); } 
    }
    for($ind = 0; $ind < $k-$out; $ind++){ 
        if($mutrt > 0){
            $mutate0 = int(double_rand((1/$mutrt))); 
            $mutate1 = int(double_rand((1/$mutrt)));
        } 
        else { $mutate0 = 0; $mutate1 = 0; }
        if($mutate0 == 42 && $mutate1 == 42){
            push(@malesB7, [ [[1, 1], [0, 0]], [[0]] ]); }
        elsif($mutate0 == 42 || $mutate1 == 42){
            push(@malesB7, [ [[1, 0], [0, 0]], [[0]] ]); }        
        else{
            push(@malesB7, [ [[0, 0], [0, 0]], [[0]] ]); } 
    }

    ### Habitat 8
    for($ind = 0; $ind < $out; $ind++){ 
        if($mutrt > 0){
            $mutate0 = int(double_rand((1/$mutrt))); 
            $mutate1 = int(double_rand((1/$mutrt)));
        } 
        else { $mutate0 = 0; $mutate1 = 0; }
        if($mutate0 == 42 && $mutate1 == 42){
            push(@femalesB8, [ [[1, 1], [0, 0]], [[0]] ]); }
        elsif($mutate0 == 42 || $mutate1 == 42){
            push(@femalesB8, [ [[1, 0], [0, 0]], [[0]] ]); }            
        else{
            push(@femalesB8, [ [[0, 0], [0, 0]], [[0]] ]); } 
    }
    for($ind = 0; $ind < $k-$out; $ind++){ 
        if($mutrt > 0){
            $mutate0 = int(double_rand((1/$mutrt))); 
            $mutate1 = int(double_rand((1/$mutrt)));
        } 
        else { $mutate0 = 0; $mutate1 = 0; }
        if($mutate0 == 42 && $mutate1 == 42){
            push(@malesB8, [ [[1, 1], [0, 0]], [[0]] ]); }
        elsif($mutate0 == 42 || $mutate1 == 42){
            push(@malesB8, [ [[1, 0], [0, 0]], [[0]] ]); }        
        else{
            push(@malesB8, [ [[0, 0], [0, 0]], [[0]] ]); } 
    }

    # ### Habitat 9
    for($ind = 0; $ind < $out; $ind++){
        if($mutrt > 0){
            $mutate0 = int(double_rand((1/$mutrt))); 
            $mutate1 = int(double_rand((1/$mutrt)));
        } 
        else { $mutate0 = 0; $mutate1 = 0; }

        if($mutate0 == 42 && $mutate1 == 42){
            push(@femalesB9, [ [[1, 1], [0, 0]], [[0]] ]); }
        elsif($mutate0 == 42 || $mutate1 == 42){
            push(@femalesB9, [ [[1, 0], [0, 0]], [[0]] ]); }            
        else{
            push(@femalesB9, [ [[0, 0], [0, 0]], [[0]] ]); } 
    }
    for($ind = 0; $ind < $k-$out; $ind++){ 
        if($mutrt > 0){
            $mutate0 = int(double_rand((1/$mutrt))); 
            $mutate1 = int(double_rand((1/$mutrt)));
        } 
        else { $mutate0 = 0; $mutate1 = 0; }
        if($mutate0 == 42 && $mutate1 == 42){
            push(@malesB9, [ [[1, 1], [0, 0]], [[0]] ]); }
        elsif($mutate0 == 42 || $mutate1 == 42){
            push(@malesB9, [ [[1, 0], [0, 0]], [[0]] ]); }        
        else{
            push(@malesB9, [ [[0, 0], [0, 0]], [[0]] ]); } 
    }

    ### Habitat 10
    for($ind = 0; $ind < $out; $ind++){
        if($mutrt > 0){
            $mutate0 = int(double_rand((1/$mutrt))); 
            $mutate1 = int(double_rand((1/$mutrt)));
        } 
        else { $mutate0 = 0; $mutate1 = 0; }
        if($mutate0 == 42 && $mutate1 == 42){
            push(@femalesB10, [ [[1, 1], [0, 0]], [[0]] ]); }
        elsif($mutate0 == 42 || $mutate1 == 42){
            push(@femalesB10, [ [[1, 0], [0, 0]], [[0]] ]); }            
        else{
            push(@femalesB10, [ [[0, 0], [0, 0]], [[0]] ]); } 
    }
    for($ind = 0; $ind < $k-$out; $ind++){ 
        if($mutrt > 0){
            $mutate0 = int(double_rand((1/$mutrt))); 
            $mutate1 = int(double_rand((1/$mutrt)));
        } 
        else { $mutate0 = 0; $mutate1 = 0; }
        if($mutate0 == 42 && $mutate1 == 42){
            push(@malesB10, [ [[1, 1], [0, 0]], [[0]] ]); }
        elsif($mutate0 == 42 || $mutate1 == 42){
            push(@malesB10, [ [[1, 0], [0, 0]], [[0]] ]); }        
        else{
            push(@malesB10, [ [[0, 0], [0, 0]], [[0]] ]); } 
    }

    
    
    ##################################################
    ############# Run through generations ############
    ##################################################

    
    # will continue the simulation until one of the species goes extinct or until endgen 
    do {


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
        
        
        
        ### Determine number of cannibal males and females in Species A and B in all Habitats
        $lenFA1c = countCannibals(@femalesA1); $lenMA1c = countCannibals(@malesA1); 
        $lenFA2c = countCannibals(@femalesA2); $lenMA2c = countCannibals(@malesA2); 
        $lenFA3c = countCannibals(@femalesA3); $lenMA3c = countCannibals(@malesA3); 
        $lenFA4c = countCannibals(@femalesA4); $lenMA4c = countCannibals(@malesA4); 
        $lenFA5c = countCannibals(@femalesA5); $lenMA5c = countCannibals(@malesA5); 
        $lenFA6c = countCannibals(@femalesA6); $lenMA6c = countCannibals(@malesA6); 
        $lenFA7c = countCannibals(@femalesA7); $lenMA7c = countCannibals(@malesA7); 
        $lenFA8c = countCannibals(@femalesA8); $lenMA8c = countCannibals(@malesA8); 
        $lenFA9c = countCannibals(@femalesA9); $lenMA9c = countCannibals(@malesA9); 
        $lenFA10c = countCannibals(@femalesA10); $lenMA10c = countCannibals(@malesA10); 

        $lenFB1c = countCannibals(@femalesB1); $lenMB1c = countCannibals(@malesB1); 
        $lenFB2c = countCannibals(@femalesB2); $lenMB2c = countCannibals(@malesB2); 
        $lenFB3c = countCannibals(@femalesB3); $lenMB3c = countCannibals(@malesB3); 
        $lenFB4c = countCannibals(@femalesB4); $lenMB4c = countCannibals(@malesB4); 
        $lenFB5c = countCannibals(@femalesB5); $lenMB5c = countCannibals(@malesB5); 
        $lenFB6c = countCannibals(@femalesB6); $lenMB6c = countCannibals(@malesB6); 
        $lenFB7c = countCannibals(@femalesB7); $lenMB7c = countCannibals(@malesB7); 
        $lenFB8c = countCannibals(@femalesB8); $lenMB8c = countCannibals(@malesB8); 
        $lenFB9c = countCannibals(@femalesB9); $lenMB9c = countCannibals(@malesB9); 
        $lenFB10c = countCannibals(@femalesB10); $lenMB10c = countCannibals(@malesB10); 
        
        ### Calculate number of non-cannibals male and females in Species A and B in all Habitats 
        $lenFA1n = $lenFA1 - $lenFA1c; $lenMA1n = $lenMA1 - $lenMA1c; 
        $lenFA2n = $lenFA2 - $lenFA2c; $lenMA2n = $lenMA2 - $lenMA2c; 
        $lenFA3n = $lenFA3 - $lenFA3c; $lenMA3n = $lenMA3 - $lenMA3c; 
        $lenFA4n = $lenFA4 - $lenFA4c; $lenMA4n = $lenMA4 - $lenMA4c; 
        $lenFA5n = $lenFA5 - $lenFA5c; $lenMA5n = $lenMA5 - $lenMA5c; 
        $lenFA6n = $lenFA6 - $lenFA6c; $lenMA6n = $lenMA6 - $lenMA6c; 
        $lenFA7n = $lenFA7 - $lenFA7c; $lenMA7n = $lenMA7 - $lenMA7c; 
        $lenFA8n = $lenFA8 - $lenFA8c; $lenMA8n = $lenMA8 - $lenMA8c; 
        $lenFA9n = $lenFA9 - $lenFA9c; $lenMA9n = $lenMA9 - $lenMA9c; 
        $lenFA10n = $lenFA10 - $lenFA10c; $lenMA10n = $lenMA10 - $lenMA10c; 

        $lenFB1n = $lenFB1 - $lenFB1c; $lenMB1n = $lenMB1 - $lenMB1c; 
        $lenFB2n = $lenFB2 - $lenFB2c; $lenMB2n = $lenMB2 - $lenMB2c; 
        $lenFB3n = $lenFB3 - $lenFB3c; $lenMB3n = $lenMB3 - $lenMB3c; 
        $lenFB4n = $lenFB4 - $lenFB4c; $lenMB4n = $lenMB4 - $lenMB4c; 
        $lenFB5n = $lenFB5 - $lenFB5c; $lenMB5n = $lenMB5 - $lenMB5c; 
        $lenFB6n = $lenFB6 - $lenFB6c; $lenMB6n = $lenMB6 - $lenMB6c; 
        $lenFB7n = $lenFB7 - $lenFB7c; $lenMB7n = $lenMB7 - $lenMB7c; 
        $lenFB8n = $lenFB8 - $lenFB8c; $lenMB8n = $lenMB8 - $lenMB8c; 
        $lenFB9n = $lenFB9 - $lenFB9c; $lenMB9n = $lenMB9 - $lenMB9c; 
        $lenFB10n = $lenFB10 - $lenFB10c; $lenMB10n = $lenMB10 - $lenMB10c; 
        
        ### Total Cannibals and Non-cannibals in Species A and B in all Habitats 
        $totalA1c = $lenFA1c + $lenMA1c; $totalA1n = $lenFA1n + $lenMA1n; 
        $totalA2c = $lenFA2c + $lenMA2c; $totalA2n = $lenFA2n + $lenMA2n; 
        $totalA3c = $lenFA3c + $lenMA3c; $totalA3n = $lenFA3n + $lenMA3n; 
        $totalA4c = $lenFA4c + $lenMA4c; $totalA4n = $lenFA4n + $lenMA4n; 
        $totalA5c = $lenFA5c + $lenMA5c; $totalA5n = $lenFA5n + $lenMA5n; 
        $totalA6c = $lenFA6c + $lenMA6c; $totalA6n = $lenFA6n + $lenMA6n; 
        $totalA7c = $lenFA7c + $lenMA7c; $totalA7n = $lenFA7n + $lenMA7n; 
        $totalA8c = $lenFA8c + $lenMA8c; $totalA8n = $lenFA8n + $lenMA8n; 
        $totalA9c = $lenFA9c + $lenMA9c; $totalA9n = $lenFA9n + $lenMA9n; 
        $totalA10c = $lenFA10c + $lenMA10c; $totalA10n = $lenFA10n + $lenMA10n; 

        $totalB1c = $lenFB1c + $lenMB1c; $totalB1n = $lenFB1n + $lenMB1n; 
        $totalB2c = $lenFB2c + $lenMB2c; $totalB2n = $lenFB2n + $lenMB2n; 
        $totalB3c = $lenFB3c + $lenMB3c; $totalB3n = $lenFB3n + $lenMB3n; 
        $totalB4c = $lenFB4c + $lenMB4c; $totalB4n = $lenFB4n + $lenMB4n; 
        $totalB5c = $lenFB5c + $lenMB5c; $totalB5n = $lenFB5n + $lenMB5n; 
        $totalB6c = $lenFB6c + $lenMB6c; $totalB6n = $lenFB6n + $lenMB6n; 
        $totalB7c = $lenFB7c + $lenMB7c; $totalB7n = $lenFB7n + $lenMB7n; 
        $totalB8c = $lenFB8c + $lenMB8c; $totalB8n = $lenFB8n + $lenMB8n; 
        $totalB9c = $lenFB9c + $lenMB9c; $totalB9n = $lenFB9n + $lenMB9n; 
        $totalB10c = $lenFB10c + $lenMB10c; $totalB10n = $lenFB10n + $lenMB10n; 
        
        $totalAc = $totalA1c + $totalA2c + $totalA3c + $totalA4c + $totalA5c + $totalA6c + $totalA7c + $totalA8c + $totalA9c + $totalA10c; 
        $totalBc = $totalB1c + $totalB2c + $totalB3c + $totalB4c + $totalB5c + $totalB6c + $totalB7c + $totalB8c + $totalB9c + $totalB10c; 

        $totalAn = $totalA - $totalAc; 
        $totalBn = $totalB - $totalBc; 
        
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
        
        makeNextGenCompetingSpecies($lenFA1, $lenMA1, $lenFB1, $lenMB1, $max_offspring1, @femalesA1, @malesA1, @femalesB1, @malesB1, @FoffspringA1, @MoffspringA1, @FoffspringB1, @MoffspringB1, @female_trackerA1, @male_trackerA1, @female_trackerB1, @male_trackerB1); 
        makeNextGenCompetingSpecies($lenFA2, $lenMA2, $lenFB2, $lenMB2, $max_offspring2, @femalesA2, @malesA2, @femalesB2, @malesB2, @FoffspringA2, @MoffspringA2, @FoffspringB2, @MoffspringB2, @female_trackerA2, @male_trackerA2, @female_trackerB2, @male_trackerB2); 
        makeNextGenCompetingSpecies($lenFA3, $lenMA3, $lenFB3, $lenMB3, $max_offspring3, @femalesA3, @malesA3, @femalesB3, @malesB3, @FoffspringA3, @MoffspringA3, @FoffspringB3, @MoffspringB3, @female_trackerA3, @male_trackerA3, @female_trackerB3, @male_trackerB3); 
        makeNextGenCompetingSpecies($lenFA4, $lenMA4, $lenFB4, $lenMB4, $max_offspring4, @femalesA4, @malesA4, @femalesB4, @malesB4, @FoffspringA4, @MoffspringA4, @FoffspringB4, @MoffspringB4, @female_trackerA4, @male_trackerA4, @female_trackerB4, @male_trackerB4); 
        makeNextGenCompetingSpecies($lenFA5, $lenMA5, $lenFB5, $lenMB5, $max_offspring5, @femalesA5, @malesA5, @femalesB5, @malesB5, @FoffspringA5, @MoffspringA5, @FoffspringB5, @MoffspringB5, @female_trackerA5, @male_trackerA5, @female_trackerB5, @male_trackerB5); 
        makeNextGenCompetingSpecies($lenFA6, $lenMA6, $lenFB6, $lenMB6, $max_offspring6, @femalesA6, @malesA6, @femalesB6, @malesB6, @FoffspringA6, @MoffspringA6, @FoffspringB6, @MoffspringB6, @female_trackerA6, @male_trackerA6, @female_trackerB6, @male_trackerB6); 
        makeNextGenCompetingSpecies($lenFA7, $lenMA7, $lenFB7, $lenMB7, $max_offspring7, @femalesA7, @malesA7, @femalesB7, @malesB7, @FoffspringA7, @MoffspringA7, @FoffspringB7, @MoffspringB7, @female_trackerA7, @male_trackerA7, @female_trackerB7, @male_trackerB7); 
        makeNextGenCompetingSpecies($lenFA8, $lenMA8, $lenFB8, $lenMB8, $max_offspring8, @femalesA8, @malesA8, @femalesB8, @malesB8, @FoffspringA8, @MoffspringA8, @FoffspringB8, @MoffspringB8, @female_trackerA8, @male_trackerA8, @female_trackerB8, @male_trackerB8); 
        makeNextGenCompetingSpecies($lenFA9, $lenMA9, $lenFB9, $lenMB9, $max_offspring9, @femalesA9, @malesA9, @femalesB9, @malesB9, @FoffspringA9, @MoffspringA9, @FoffspringB9, @MoffspringB9, @female_trackerA9, @male_trackerA9, @female_trackerB9, @male_trackerB9); 
        makeNextGenCompetingSpecies($lenFA10, $lenMA10, $lenFB10, $lenMB10, $max_offspring10, @femalesA10, @malesA10, @femalesB10, @malesB10, @FoffspringA10, @MoffspringA10, @FoffspringB10, @MoffspringB10, @female_trackerA10, @male_trackerA10, @female_trackerB10, @male_trackerB10); 
               
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
            migration(1, 'A', @femalesA1, @female_trackerA1, @male_trackerA1, @FmigrantsA1, @FmigrantsA2, @FmigrantsA3, @FmigrantsA4, @FmigrantsA5, @FmigrantsA6, @FmigrantsA7, @FmigrantsA8, @FmigrantsA9, @FmigrantsA10, 
                                @MmigrantsA1, @MmigrantsA2, @MmigrantsA3, @MmigrantsA4, @MmigrantsA5, @MmigrantsA6, @MmigrantsA7, @MmigrantsA8, @MmigrantsA9, @MmigrantsA10, 
                                @FoffspringA1, @MoffspringA1);
                                
            migration(2, 'A', @femalesA2, @female_trackerA2, @male_trackerA2, @FmigrantsA1, @FmigrantsA2, @FmigrantsA3, @FmigrantsA4, @FmigrantsA5, @FmigrantsA6, @FmigrantsA7, @FmigrantsA8, @FmigrantsA9, @FmigrantsA10, 
                                @MmigrantsA1, @MmigrantsA2, @MmigrantsA3, @MmigrantsA4, @MmigrantsA5, @MmigrantsA6, @MmigrantsA7, @MmigrantsA8, @MmigrantsA9, @MmigrantsA10, 
                                @FoffspringA2, @MoffspringA2);
                                
            migration(3, 'A', @femalesA3, @female_trackerA3, @male_trackerA3, @FmigrantsA1, @FmigrantsA2, @FmigrantsA3, @FmigrantsA4, @FmigrantsA5, @FmigrantsA6, @FmigrantsA7, @FmigrantsA8, @FmigrantsA9, @FmigrantsA10, 
                                @MmigrantsA1, @MmigrantsA2, @MmigrantsA3, @MmigrantsA4, @MmigrantsA5, @MmigrantsA6, @MmigrantsA7, @MmigrantsA8, @MmigrantsA9, @MmigrantsA10, 
                                @FoffspringA3, @MoffspringA3);
                                
            migration(4, 'A', @femalesA4, @female_trackerA4, @male_trackerA4, @FmigrantsA1, @FmigrantsA2, @FmigrantsA3, @FmigrantsA4, @FmigrantsA5, @FmigrantsA6, @FmigrantsA7, @FmigrantsA8, @FmigrantsA9, @FmigrantsA10, 
                                @MmigrantsA1, @MmigrantsA2, @MmigrantsA3, @MmigrantsA4, @MmigrantsA5, @MmigrantsA6, @MmigrantsA7, @MmigrantsA8, @MmigrantsA9, @MmigrantsA10, 
                                @FoffspringA4, @MoffspringA4);
                                
            migration(5, 'A', @femalesA5, @female_trackerA5, @male_trackerA5, @FmigrantsA1, @FmigrantsA2, @FmigrantsA3, @FmigrantsA4, @FmigrantsA5, @FmigrantsA6, @FmigrantsA7, @FmigrantsA8, @FmigrantsA9, @FmigrantsA10, 
                                @MmigrantsA1, @MmigrantsA2, @MmigrantsA3, @MmigrantsA4, @MmigrantsA5, @MmigrantsA6, @MmigrantsA7, @MmigrantsA8, @MmigrantsA9, @MmigrantsA10, 
                                @FoffspringA5, @MoffspringA5);
                                
            migration(6, 'A', @femalesA6, @female_trackerA6, @male_trackerA6, @FmigrantsA1, @FmigrantsA2, @FmigrantsA3, @FmigrantsA4, @FmigrantsA5, @FmigrantsA6, @FmigrantsA7, @FmigrantsA8, @FmigrantsA9, @FmigrantsA10, 
                                @MmigrantsA1, @MmigrantsA2, @MmigrantsA3, @MmigrantsA4, @MmigrantsA5, @MmigrantsA6, @MmigrantsA7, @MmigrantsA8, @MmigrantsA9, @MmigrantsA10, 
                                @FoffspringA6, @MoffspringA6);
                                
            migration(7, 'A', @femalesA7, @female_trackerA7, @male_trackerA7, @FmigrantsA1, @FmigrantsA2, @FmigrantsA3, @FmigrantsA4, @FmigrantsA5, @FmigrantsA6, @FmigrantsA7, @FmigrantsA8, @FmigrantsA9, @FmigrantsA10, 
                                @MmigrantsA1, @MmigrantsA2, @MmigrantsA3, @MmigrantsA4, @MmigrantsA5, @MmigrantsA6, @MmigrantsA7, @MmigrantsA8, @MmigrantsA9, @MmigrantsA10, 
                                @FoffspringA7, @MoffspringA7);
                                
            migration(8, 'A', @femalesA8, @female_trackerA8, @male_trackerA8, @FmigrantsA1, @FmigrantsA2, @FmigrantsA3, @FmigrantsA4, @FmigrantsA5, @FmigrantsA6, @FmigrantsA7, @FmigrantsA8, @FmigrantsA9, @FmigrantsA10, 
                                @MmigrantsA1, @MmigrantsA2, @MmigrantsA3, @MmigrantsA4, @MmigrantsA5, @MmigrantsA6, @MmigrantsA7, @MmigrantsA8, @MmigrantsA9, @MmigrantsA10, 
                                @FoffspringA8, @MoffspringA8);
                                
            migration(9, 'A', @femalesA9, @female_trackerA9, @male_trackerA9, @FmigrantsA1, @FmigrantsA2, @FmigrantsA3, @FmigrantsA4, @FmigrantsA5, @FmigrantsA6, @FmigrantsA7, @FmigrantsA8, @FmigrantsA9, @FmigrantsA10, 
                                @MmigrantsA1, @MmigrantsA2, @MmigrantsA3, @MmigrantsA4, @MmigrantsA5, @MmigrantsA6, @MmigrantsA7, @MmigrantsA8, @MmigrantsA9, @MmigrantsA10, 
                                @FoffspringA9, @MoffspringA9);
                                
            migration(10, 'A', @femalesA10, @female_trackerA10, @male_trackerA10, @FmigrantsA1, @FmigrantsA2, @FmigrantsA3, @FmigrantsA4, @FmigrantsA5, @FmigrantsA6, @FmigrantsA7, @FmigrantsA8, @FmigrantsA9, @FmigrantsA10, 
                                @MmigrantsA1, @MmigrantsA2, @MmigrantsA3, @MmigrantsA4, @MmigrantsA5, @MmigrantsA6, @MmigrantsA7, @MmigrantsA8, @MmigrantsA9, @MmigrantsA10, 
                                @FoffspringA10, @MoffspringA10);
                                
           
            migration(1, 'B', @femalesB1, @female_trackerB1, @male_trackerB1, @FmigrantsB1, @FmigrantsB2, @FmigrantsB3, @FmigrantsB4, @FmigrantsB5, @FmigrantsB6, @FmigrantsB7, @FmigrantsB8, @FmigrantsB9, @FmigrantsB10, 
                                @MmigrantsB1, @MmigrantsB2, @MmigrantsB3, @MmigrantsB4, @MmigrantsB5, @MmigrantsB6, @MmigrantsB7, @MmigrantsB8, @MmigrantsB9, @MmigrantsB10, 
                                @FoffspringB1, @MoffspringB1);
            
            migration(2, 'B', @femalesB2, @female_trackerB2, @male_trackerB2, @FmigrantsB1, @FmigrantsB2, @FmigrantsB3, @FmigrantsB4, @FmigrantsB5, @FmigrantsB6, @FmigrantsB7, @FmigrantsB8, @FmigrantsB9, @FmigrantsB10, 
                                @MmigrantsB1, @MmigrantsB2, @MmigrantsB3, @MmigrantsB4, @MmigrantsB5, @MmigrantsB6, @MmigrantsB7, @MmigrantsB8, @MmigrantsB9, @MmigrantsB10, 
                                @FoffspringB2, @MoffspringB2);
            
            migration(3, 'B', @femalesB3, @female_trackerB3, @male_trackerB3, @FmigrantsB1, @FmigrantsB2, @FmigrantsB3, @FmigrantsB4, @FmigrantsB5, @FmigrantsB6, @FmigrantsB7, @FmigrantsB8, @FmigrantsB9, @FmigrantsB10, 
                                @MmigrantsB1, @MmigrantsB2, @MmigrantsB3, @MmigrantsB4, @MmigrantsB5, @MmigrantsB6, @MmigrantsB7, @MmigrantsB8, @MmigrantsB9, @MmigrantsB10, 
                                @FoffspringB3, @MoffspringB3);
            
            migration(4, 'B', @femalesB4, @female_trackerB4, @male_trackerB4, @FmigrantsB1, @FmigrantsB2, @FmigrantsB3, @FmigrantsB4, @FmigrantsB5, @FmigrantsB6, @FmigrantsB7, @FmigrantsB8, @FmigrantsB9, @FmigrantsB10, 
                                @MmigrantsB1, @MmigrantsB2, @MmigrantsB3, @MmigrantsB4, @MmigrantsB5, @MmigrantsB6, @MmigrantsB7, @MmigrantsB8, @MmigrantsB9, @MmigrantsB10, 
                                @FoffspringB4, @MoffspringB4);
            
            migration(5, 'B', @femalesB5, @female_trackerB5, @male_trackerB5, @FmigrantsB1, @FmigrantsB2, @FmigrantsB3, @FmigrantsB4, @FmigrantsB5, @FmigrantsB6, @FmigrantsB7, @FmigrantsB8, @FmigrantsB9, @FmigrantsB10, 
                                @MmigrantsB1, @MmigrantsB2, @MmigrantsB3, @MmigrantsB4, @MmigrantsB5, @MmigrantsB6, @MmigrantsB7, @MmigrantsB8, @MmigrantsB9, @MmigrantsB10, 
                                @FoffspringB5, @MoffspringB5);
            
            migration(6, 'B', @femalesB6, @female_trackerB6, @male_trackerB6, @FmigrantsB1, @FmigrantsB2, @FmigrantsB3, @FmigrantsB4, @FmigrantsB5, @FmigrantsB6, @FmigrantsB7, @FmigrantsB8, @FmigrantsB9, @FmigrantsB10, 
                                @MmigrantsB1, @MmigrantsB2, @MmigrantsB3, @MmigrantsB4, @MmigrantsB5, @MmigrantsB6, @MmigrantsB7, @MmigrantsB8, @MmigrantsB9, @MmigrantsB10, 
                                @FoffspringB6, @MoffspringB6);
            
            migration(7, 'B', @femalesB7, @female_trackerB7, @male_trackerB7, @FmigrantsB1, @FmigrantsB2, @FmigrantsB3, @FmigrantsB4, @FmigrantsB5, @FmigrantsB6, @FmigrantsB7, @FmigrantsB8, @FmigrantsB9, @FmigrantsB10, 
                                @MmigrantsB1, @MmigrantsB2, @MmigrantsB3, @MmigrantsB4, @MmigrantsB5, @MmigrantsB6, @MmigrantsB7, @MmigrantsB8, @MmigrantsB9, @MmigrantsB10, 
                                @FoffspringB7, @MoffspringB7);
            
            migration(8, 'B', @femalesB8, @female_trackerB8, @male_trackerB8, @FmigrantsB1, @FmigrantsB2, @FmigrantsB3, @FmigrantsB4, @FmigrantsB5, @FmigrantsB6, @FmigrantsB7, @FmigrantsB8, @FmigrantsB9, @FmigrantsB10, 
                                @MmigrantsB1, @MmigrantsB2, @MmigrantsB3, @MmigrantsB4, @MmigrantsB5, @MmigrantsB6, @MmigrantsB7, @MmigrantsB8, @MmigrantsB9, @MmigrantsB10, 
                                @FoffspringB8, @MoffspringB8);
            
            migration(9, 'B', @femalesB9, @female_trackerB9, @male_trackerB9, @FmigrantsB1, @FmigrantsB2, @FmigrantsB3, @FmigrantsB4, @FmigrantsB5, @FmigrantsB6, @FmigrantsB7, @FmigrantsB8, @FmigrantsB9, @FmigrantsB10, 
                                @MmigrantsB1, @MmigrantsB2, @MmigrantsB3, @MmigrantsB4, @MmigrantsB5, @MmigrantsB6, @MmigrantsB7, @MmigrantsB8, @MmigrantsB9, @MmigrantsB10, 
                                @FoffspringB9, @MoffspringB9);
            
            migration(10, 'B', @femalesB10, @female_trackerB10, @male_trackerB10, @FmigrantsB1, @FmigrantsB2, @FmigrantsB3, @FmigrantsB4, @FmigrantsB5, @FmigrantsB6, @FmigrantsB7, @FmigrantsB8, @FmigrantsB9, @FmigrantsB10, 
                                @MmigrantsB1, @MmigrantsB2, @MmigrantsB3, @MmigrantsB4, @MmigrantsB5, @MmigrantsB6, @MmigrantsB7, @MmigrantsB8, @MmigrantsB9, @MmigrantsB10, 
                                @FoffspringB10, @MoffspringB10);
                                
                                              
        
        
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
  
        
        
        print(FH2 "$totalA1c,$totalA1n,$totalA2c,$totalA2n,$totalA3c,$totalA3n,$totalA4c,$totalA4n,$totalA5c,$totalA5n,$totalA6c,$totalA6n,$totalA7c,$totalA7n,$totalA8c,$totalA8n,$totalA9c,$totalA9n,$totalA10c,$totalA10n,$totalB1c,$totalB1n,$totalB2c,$totalB2n,$totalB3c,$totalB3n,$totalB4c,$totalB4n,$totalB5c,$totalB5n,$totalB6c,$totalB6n,$totalB7c,$totalB7n,$totalB8c,$totalB8n,$totalB9c,$totalB9n,$totalB10c,$totalB10n,$chance_cannibalism_cost,$chance_cannibalism_offspring,$gen,$simno\n"); 
        print(FH3 "$totalA1c,$totalA1n,$totalA2c,$totalA2n,$totalA3c,$totalA3n,$totalA4c,$totalA4n,$totalA5c,$totalA5n,$totalA6c,$totalA6n,$totalA7c,$totalA7n,$totalA8c,$totalA8n,$totalA9c,$totalA9n,$totalA10c,$totalA10n,$totalB1c,$totalB1n,$totalB2c,$totalB2n,$totalB3c,$totalB3n,$totalB4c,$totalB4n,$totalB5c,$totalB5n,$totalB6c,$totalB6n,$totalB7c,$totalB7n,$totalB8c,$totalB8n,$totalB9c,$totalB9n,$totalB10c,$totalB10n,$gen,$simno\n"); 

if ($gen==$endgen){
        print(FH4 "$gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,A,1,no,$totalA1n,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,A,1,yes,$totalA1c,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,A,2,no,$totalA2n,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,A,2,yes,$totalA2c,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,A,3,no,$totalA3n,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,A,3,yes,$totalA3c,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,A,4,no,$totalA4n,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,A,4,yes,$totalA4c,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,A,5,no,$totalA5n,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,A,5,yes,$totalA5c,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,A,6,no,$totalA6n,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,A,6,yes,$totalA6c,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,A,7,no,$totalA7n,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,A,7,yes,$totalA7c,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,A,8,no,$totalA8n,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,A,8,yes,$totalA8c,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,A,9,no,$totalA9n,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,A,9,yes,$totalA9c,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,A,10,no,$totalA10n,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,A,10,yes,$totalA10c,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,B,1,no,space,$totalB1n,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,B,1,yes,space,$totalB1c,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,B,2,no,space,$totalB2n,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,B,2,yes,space,$totalB2c,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,B,3,no,space,$totalB3n,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,B,3,yes,space,$totalB3c,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,B,4,no,space,$totalB4n,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,B,4,yes,space,$totalB4c,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,B,5,no,space,$totalB5n,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,B,5,yes,space,$totalB5c,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,B,6,no,space,$totalB6n,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,B,6,yes,space,$totalB6c,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,B,7,no,space,$totalB7n,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,B,7,yes,space,$totalB7c,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,B,8,no,space,$totalB8n,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,B,8,yes,space,$totalB8c,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,B,9,no,space,$totalB9n,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,B,9,yes,space,$totalB9c,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,B,10,no,space,$totalB10n,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,B,10,yes,space,$totalB10c\n");
		print(FH5"$simno,$propA1c,$totalA1,$propB1c,$totalB1\n,$propA2c,$totalA2,$propB2c,$totalB2\n,$propA3c,$totalA3,$propB3c,$totalB3\n,$propA4c,$totalA4,$propB4c,$totalB4\n,$propA5c,$totalA5,$propB5c,$totalB5\n,$propA6c,$totalA6,$propB6c,$totalB6\n,$propA7c,$totalA7,$propB7c,$totalB7\n,$propA8c,$totalA8,$propB8c,$totalB8\n,$propA9c,$totalA9,$propB9c,$totalB9\n,$propA10c,$totalA10,$propB10c,$totalB10\n");
    
	print(FH4d "$gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,A,1,no,$totalA1n,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,A,1,yes,$totalA1c,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,A,2,no,$totalA2n,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,A,2,yes,$totalA2c,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,A,3,no,$totalA3n,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,A,3,yes,$totalA3c,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,A,4,no,$totalA4n,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,A,4,yes,$totalA4c,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,A,5,no,$totalA5n,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,A,5,yes,$totalA5c,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,A,6,no,$totalA6n,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,A,6,yes,$totalA6c,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,A,7,no,$totalA7n,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,A,7,yes,$totalA7c,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,A,8,no,$totalA8n,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,A,8,yes,$totalA8c,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,A,9,no,$totalA9n,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,A,9,yes,$totalA9c,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,A,10,no,$totalA10n,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,A,10,yes,$totalA10c,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,B,1,no,space,$totalB1n,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,B,1,yes,space,$totalB1c,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,B,2,no,space,$totalB2n,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,B,2,yes,space,$totalB2c,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,B,3,no,space,$totalB3n,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,B,3,yes,space,$totalB3c,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,B,4,no,space,$totalB4n,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,B,4,yes,space,$totalB4c,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,B,5,no,space,$totalB5n,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,B,5,yes,space,$totalB5c,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,B,6,no,space,$totalB6n,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,B,6,yes,space,$totalB6c,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,B,7,no,space,$totalB7n,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,B,7,yes,space,$totalB7c,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,B,8,no,space,$totalB8n,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,B,8,yes,space,$totalB8c,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,B,9,no,space,$totalB9n,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,B,9,yes,space,$totalB9c,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,B,10,no,space,$totalB10n,
        $gen,$simno,$chance_cannibalism_cost,$chance_cannibalism_offspring,B,10,yes,space,$totalB10c\n");
		print(FH5d"$simno,$propA1c,$totalA1,$propB1c,$totalB1\n,$propA2c,$totalA2,$propB2c,$totalB2\n,$propA3c,$totalA3,$propB3c,$totalB3\n,$propA4c,$totalA4,$propB4c,$totalB4\n,$propA5c,$totalA5,$propB5c,$totalB5\n,$propA6c,$totalA6,$propB6c,$totalB6\n,$propA7c,$totalA7,$propB7c,$totalB7\n,$propA8c,$totalA8,$propB8c,$totalB8\n,$propA9c,$totalA9,$propB9c,$totalB9\n,$propA10c,$totalA10,$propB10c,$totalB10\n");

	}
}    

    ### The below "until" statement (other half of the "do" statement earlier) keeps the simulation going until the ending conditions are met. 
    ### Here, it's when endgen has been reached 
    until($gen >= $endgen); 

    
   $finalgen = $gen; #the number of generations the simulation ran for 
    
   if($totalA == 0) { $extinctionA = 1; $total_extinctionsA++; }
   else { $extinctionA = 0; } 

   if($totalB == 0) { $extinctionB = 1; $total_extinctionsB++; }
   else { $extinctionB = 0; } 



    ### Print results to per-simulation file
    print(FH1 "$total,$totalA,$totalB,$totalA1,$totalA2,$totalA3,$totalA4,$totalA5,$totalA6,$totalA7,$totalA8,$totalA9,$totalA10,$totalB1,$totalB2,$totalB3,$totalB4,$totalB5,$totalB6,$totalB7,$totalB8,$totalB9,$totalB10,$propATotal,$propBTotal,$finalgen,$simno\n");
print(FH1d "$total,$totalA,$totalB,$totalA1,$totalA2,$totalA3,$totalA4,$totalA5,$totalA6,$totalA7,$totalA8,$totalA9,$totalA10,$totalB1,$totalB2,$totalB3,$totalB4,$totalB5,$totalB6,$totalB7,$totalB8,$totalB9,$totalB10,$propA1c,$propA2c,$propA3c,$propA4c,$propA5c,$propA6c,$propA7c,$propA8c,$propA9c,$propA10c,$propB1c,$propB2c,$propB3c,$propB4c,$propB5c,$propB6c,$propB7c,$propB8c,$propB9c,$propB10c,$propATotal,$propBTotal,$finalgen,$simno\n");

    
    # accumulate results for summary file 
    $avg_total = $avg_total + $total; 
    $avg_total_A = $avg_total_A + $totalA; 
    $avg_total_B = $avg_total_B + $totalB; 
    $avg_prop_A = $avg_prop_A + $propATotal; 
    $avg_prop_B = $avg_prop_B + $propBTotal; 
    $avg_endgen = $avg_endgen + $finalgen; 
	$avg_propAc=$avg_propAc + $propAcTotal;
	$avg_propBc=$avg_propBc + $propBcTotal;	
    
    $avgA1 = $avgA1 + $totalA1; $avgB1 = $avgB1 + $totalB1; 
    $avgA2 = $avgA2 + $totalA2; $avgB2 = $avgB2 + $totalB2; 
    $avgA3 = $avgA3 + $totalA3; $avgB3 = $avgB3 + $totalB3; 
    $avgA4 = $avgA4 + $totalA4; $avgB4 = $avgB4 + $totalB4; 
    $avgA5 = $avgA5 + $totalA5; $avgB5 = $avgB5 + $totalB5; 
    $avgA6 = $avgA6 + $totalA6; $avgB6 = $avgB6 + $totalB6; 
    $avgA7 = $avgA7 + $totalA7; $avgB7 = $avgB7 + $totalB7; 
    $avgA8 = $avgA8 + $totalA8; $avgB8 = $avgB8 + $totalB8; 
    $avgA9 = $avgA9 + $totalA9; $avgB9 = $avgB9 + $totalB9; 
    $avgA10 = $avgA10 + $totalA10; $avgB10 = $avgB10 + $totalB10; 
    
        if ($propA1c>0){$cantalA++};
		if ($propA2c>0){$cantalA++};
		if ($propA3c>0){$cantalA++};
		if ($propA4c>0){$cantalA++};
		if ($propA5c>0){$cantalA++};
		if ($propA6c>0){$cantalA++};
		if ($propA7c>0){$cantalA++};
		if ($propA8c>0){$cantalA++};
		if ($propA9c>0){$cantalA++};
		if ($propA10c>0){$cantalA++};
		if ($propB1c>0){$cantalB++};
		if ($propB2c>0){$cantalB++};
		if ($propB3c>0){$cantalB++};
		if ($propB4c>0){$cantalB++};
		if ($propB5c>0){$cantalB++};
		if ($propB6c>0){$cantalB++};
		if ($propB7c>0){$cantalB++};
		if ($propB8c>0){$cantalB++};
		if ($propB9c>0){$cantalB++};
		if ($propB10c>0){$cantalB++};

        if ($propA1c>0.5){$cantal50A++};
		if ($propA2c>0.5){$cantal50A++};
		if ($propA3c>0.5){$cantal50A++};
		if ($propA4c>0.5){$cantal50A++};
		if ($propA5c>0.5){$cantal50A++};
		if ($propA6c>0.5){$cantal50A++};
		if ($propA7c>0.5){$cantal50A++};
		if ($propA8c>0.5){$cantal50A++};
		if ($propA9c>0.5){$cantal50A++};
		if ($propA10c>0.5){$cantal50A++};
		if ($propB1c>0.5){$cantal50B++};
		if ($propB2c>0.5){$cantal50B++};
		if ($propB3c>0.5){$cantal50B++};
		if ($propB4c>0.5){$cantal50B++};
		if ($propB5c>0.5){$cantal50B++};
		if ($propB6c>0.5){$cantal50B++};
		if ($propB7c>0.5){$cantal50B++};
		if ($propB8c>0.5){$cantal50B++};
		if ($propB9c>0.5){$cantal50B++};
		if ($propB10c>0.5){$cantal50B++};
		
        if ($totalA1==0){$extalA++,$exorA++};
		if ($totalA2==0){$extalA++,$exorA++};
		if ($totalA3==0){$extalA++,$exorA++};
		if ($totalA4==0){$extalA++,$exorA++};
		if ($totalA5==0){$extalA++,$exorA++};
		if ($totalA6==0){$extalA++};
		if ($totalA7==0){$extalA++};
		if ($totalA8==0){$extalA++};
		if ($totalA9==0){$extalA++};
		if ($totalA10==0){$extalA++};
		if ($totalB1==0){$extalB++};
		if ($totalB2==0){$extalB++};
		if ($totalB3==0){$extalB++};
		if ($totalB4==0){$extalB++};
		if ($totalB5==0){$extalB++};
		if ($totalB6==0){$extalB++,$exorB++};
		if ($totalB7==0){$extalB++,$exorB++};
		if ($totalB8==0){$extalB++,$exorB++};
		if ($totalB9==0){$extalB++,$exorB++};
		if ($totalB10==0){$extalB++,$exorB++};	
	
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
$avg_propAc = $avg_propAc/$runs; 
$avg_propBc = $avg_propBc/$runs; 


$avgA1 = $avgA1/$runs; $avgB1 = $avgB1/$runs; 
$avgA2 = $avgA2/$runs; $avgB2 = $avgB2/$runs; 
$avgA3 = $avgA3/$runs; $avgB3 = $avgB3/$runs; 
$avgA4 = $avgA4/$runs; $avgB4 = $avgB4/$runs; 
$avgA5 = $avgA5/$runs; $avgB5 = $avgB5/$runs; 
$avgA6 = $avgA6/$runs; $avgB6 = $avgB6/$runs; 
$avgA7 = $avgA7/$runs; $avgB7 = $avgB7/$runs; 
$avgA8 = $avgA8/$runs; $avgB8 = $avgB8/$runs; 
$avgA9 = $avgA9/$runs; $avgB9 = $avgB9/$runs; 
$avgA10 = $avgA10/$runs; $avgB10 = $avgB10/$runs; 

$cantal=$cantalA+$cantalB;
$cantal50=$cantal50A+$cantal50B;
$extaltot=$extalA+$extalB; 
$empty= $extaltot-10*$runs;
### Print cumulative results to summary file
print(FH "$avg_total,$avg_total_A,$avg_total_B,$avgA1,$avgA2,$avgA3,$avgA4,$avgA5,$avgA6,$avgA7,$avgA8,$avgA9,$avgA10,$avgB1,$avgB2,$avgB3,$avgB4,$avgB5,$avgB6,$avgB7,$avgB8,$avgB9,$avgB10,$avg_prop_A,$avg_prop_B,$avg_endgen,$runs\n");

print(FHd "$avg_total,$avg_total_A,$avg_total_B,$avgA1,$avgA2,$avgA3,$avgA4,$avgA5,$avgA6,$avgA7,$avgA8,$avgA9,$avgA10,$avgB1,$avgB2,$avgB3,$avgB4,$avgB5,$avgB6,$avgB7,$avgB8,$avgB9,$avgB10,$avg_prop_A,$avg_prop_B,$avg_endgen,$runs,$avg_propAc,$avg_propBc,$cantal,$empty\n");

# End that particular set of parameters; move on to the next set of parameters (if looping, else simulation is over) 

}   # Closing bracket for file looping 
}

# } # Closing bracket for mutation rate loop
# } # For the migration rate loop
# } # For cancost 
# } # For num_can_offspring
# } # For carrying capacity 
# } # For dead 

print (FH6 "$test,$bene,$cantal,$cantal50,$empty,$extalA,$extalB,$total_extinctionsA,$total_extinctionsB,$extorA,$extorB\n");
}}
