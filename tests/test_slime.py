
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

    def test_simulate_recent_history(self):
        slime.simulate_recent_history("ex_script.slim")

    def test_slime_run(self):
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
                growth_rate = 0),
            # ASW - needed as a 'dummy'
            msprime.PopulationConfiguration(
                sample_size = 0,                                                              
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

        mysim = slime.AdmixtureSimulation("ex_script.slim", slim_out="ex_slim.trees",
                    populations_to_sample_from = [2], sample_sizes = [100],
                    neutral_mutation_rate = mu, ancient_recombination_rate = rho,
                    ancient_population_configurations = population_configurations,
                    ancient_demographic_events = demographic_events, 
                    out_file = "ex_out.trees")
        ts = mysim.go()

