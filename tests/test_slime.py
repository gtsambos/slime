
import slime
import unittest
import msprime, pyslim

class TestSlime(unittest.TestCase):

    # Example supplied SLiM script
    slim_script = '''
// set up a simple neutral simulation
initialize() {
    
    // Set constants.
    defineConstant("chromLength", 9999);
    defineConstant("p0Size", 100);
    defineConstant("p1Size", 100);
    defineConstant("p2Size", 100);
    defineConstant("p0Prop", 1/100);
    defineConstant("p1Prop", 1 - p0Prop);

    initializeSLiMModelType("WF");
    initializeSLiMOptions(keepPedigrees = T); // Keep pedigree information
    initializeTreeSeq();
    initializeMutationRate(1e-3);
    
    // m1 mutation type: slightly positive
    initializeMutationType("m1", 0.5, "f", 0.1);
    
    // g1 genomic element type: uses m1 for all mutations
    initializeGenomicElementType("g1", m1, 1.0);
    
    // uniform chromosome of length 100 kb
    initializeGenomicElement(g1, 0, chromLength);
    
    initializeRecombinationRate(1e-8);
    
}

1 early(){
    sim.addSubpop("p0", p0Size);
    sim.addSubpop("p1", p1Size);
}

1 late() { 
    
    // Add admixing population.
    sim.addSubpop("p2", p2Size);
    p2.setMigrationRates(c(p0,p1),c(p0Prop, p1Prop));
    
}

2 late() {
    p2.setMigrationRates(c(p0,p1),c(0.0,0.0));
}

// save the tree output after 25 generations.
25 late() {
    sim.treeSeqOutput("examples/ex-recent-history.trees");
    sim.simulationFinished();
}
'''

    slim_file = open("ex_script.slim", "w")
    slim_file.writelines(slim_script)
    slim_file.close()

    scr = slime.RecentHistory(final_gen = 10, chrom_length = 10)
    scr.initialize_recombination(0.2)
    config = msprime.PopulationConfiguration(0, 100, growth_rate = 0)
    config1 = msprime.PopulationConfiguration(0, 100, growth_rate = 0.05)
    scr.add_reference_population(config, 'pop0')
    scr.add_reference_population(config, 'pop1')
    scr.add_admixed_population(config1, 'adm', proportions = [0.2, 0.8])
    change1 = msprime.PopulationParametersChange(3, growth_rate =  .5, population_id = 0)
    change2 = msprime.PopulationParametersChange(6, growth_rate =  .6, population_id = 0)
    change3 = msprime.PopulationParametersChange(2, growth_rate =  .3, population_id = 1)
    change4 = msprime.PopulationParametersChange(8, growth_rate =  0, population_id = 1, initial_size = 100)
    scr.add_demographic_events([change3, change1, change2, change4])
    scr.save_script("ex_script2.slim")

    def test_simulate_recent_history(self):
        slime.simulate_recent_history("ex_script.slim")
        slime.simulate_recent_history("ex_script2.slim")

    def test_slime_run_saved_script(self):
        rho = 1e-8
        mu = 1e-8
        ######################
        # POPULATION HISTORY #
        ######################
        N0 = 1e4 # population size
        divergenceGen = 20000

        # Population IDs
        population_configurations = [
            # CEU
            msprime.PopulationConfiguration(
                initial_size = N0,
                growth_rate = 0),
            # YRI
            msprime.PopulationConfiguration(                                
                initial_size = N0,                                
                growth_rate = 0)
        ]

        demographic_events = [
            msprime.MassMigration(
                time=divergenceGen,
                source = 1,
                destination = 0,
                proportion = 1.0),
            msprime.MigrationRateChange(time = divergenceGen,
                rate = 0)
        ]

        # First, with the created script.
        mysim = slime.AdmixtureSimulation(slim_script="ex_script.slim", 
                    slim_out="examples/ex-recent-history.trees",
                    neutral_mutation_rate = mu, ancient_recombination_rate = rho,
                    ancient_population_configurations = population_configurations,
                    ancient_demographic_events = demographic_events,
                    out_file = "examples/ex_out.trees",
                    populations_to_sample_from = [2], sample_sizes = [10], 
                    )
        ts = mysim.go()

        # Next, with the generated script.
        mysim2 = slime.AdmixtureSimulation(slim_script="ex_script2.slim", 
            slim_out="examples/ex-recent-history.trees",
            neutral_mutation_rate = mu, ancient_recombination_rate = rho,
            ancient_population_configurations = population_configurations,
            ancient_demographic_events = demographic_events,
            out_file = "examples/ex_out.trees",
            populations_to_sample_from = [2], sample_sizes = [10], 
            )
        ts2 = mysim2.go()


