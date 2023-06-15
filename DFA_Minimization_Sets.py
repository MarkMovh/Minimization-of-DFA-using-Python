import copy

# Notes:
# 1. Add Myhill-Nerode method (At the moment only conducts Equivalence Method)
# x. Adjust how the finite automata is input (ideally image) 
#       -> First be able to detect states and transitions
#       -> Be able to identify state names and transition symbols
#       -> Connect the states and transitions to build the dfa dictionary


class InvalidDeterministicFiniteAutomata(Exception):
    "Invalid Finite Automata Input Detected"
    pass

class DFA_Scanner:
    # Method ascertains that the input fa is a valid deterministic autom
    def check_dfa(self, fa, symbols):
        # Check that the length of transitions is no more than the number of symbols in the alphabet
        # If this does not hold true, raise an exception to notify user
        try:
            for transitions in fa:
                if len(fa[transitions]) == 0:
                    pass
                
                if len(fa[transitions]) > len(symbols):
                    raise InvalidDeterministicFiniteAutomata

                else:
                    pass

            return True

        except:
            raise Exception('''
        Error: The following automata you have input is not deterministic.

        An Automata can only be minimized if it is in deterministic form.
        Convert the input automata into DFA first.
                    ''')



# DFA minimizer method based on equivalence
class DFA_Equivalence_Minimizer(DFA_Scanner):
    # Initialize the class by getting the minimized set and dfa
    def __init__(self, dfa, alphabet):
        self.dfa = dfa
        self.alphabet = alphabet

        self.check_dfa()

        get_first_split = self.initial_split(self.dfa, self.alphabet)
        self.minimized_set = self.minimize_dfa(self.dfa, get_first_split, self.alphabet)
        self.minimized_dfa = self.get_transition_table(self.dfa, self.minimized_set, self.alphabet)

    # Create get methods for the initialized information
    def get_minimized_set(self):
        return self.minimized_set
    
    def get_minimized_dfa(self):
        return self.minimized_dfa

    # Using inheritance, a check method can be created to ensure the FA is DFA
    def check_dfa(self):
        return DFA_Scanner.check_dfa(self, self.dfa, self.alphabet)

    # Get all the states and separate them into two sets:
    # - One containing the final states
    # - One containing all other states
    def initial_split(self, dfa, alphabet):
        initial_set = [[], []]

        # The split is done by checking the initial of the state for the symbol *
        for state in dfa.keys():
            if state[0] == "*":
                initial_set[0].append(state)
            else:
                initial_set[1].append(state)

        return initial_set 

    # The transitions of each state are checked, and summarised in a dictionary
    def check_transitions(self, dfa, full_set, inner_set, alphabet):
        temp_tran_loc = {}

        # Get the state of the current subset
        for state in inner_set:
            symbol_list = []

            # For every alphabet symbol, collect the transitions and then check in what list they are
            for symbol in range(len(alphabet)):
                current_transition = dfa[state][alphabet[symbol]]
                for index, split in enumerate(full_set):
                    if current_transition in split:
                        # Save this list location in the symbol list
                        symbol_list.append(index)

            # For every state save their list location
            temp_tran_loc[state] = symbol_list

        return temp_tran_loc

    # The following method determines which states should be moved from their list.
    def identify_outliers(self, transitions, previous_set):
        count_dict = {}

        # Iterate over every state and its list locations, saving them to the dictionary
        for state in transitions.values():
            count_dict[tuple(state)] = count_dict.get(tuple(state), 0) + 1

        # Extract the list with the highest count to identify the majority
        majority_list = max(count_dict, key=count_dict.get)

        # Get the states that are in the minority list
        outlier_states = [key for key, value in transitions.items() if tuple(value) != majority_list]

        # Using the outlier states, check the input set for these states and then remove them.
        # The remove set collects these, where it is then appended to the end of the input set as a new set of states.
        updated_set = copy.deepcopy(previous_set)
        if len(outlier_states) > 0:
            remove_set = []
            for index, state_set in enumerate(updated_set):
                for outlier in outlier_states:
                    if outlier in state_set:
                        state_set.remove(outlier)
                        remove_set.append(outlier)
            
            updated_set.append(remove_set)

            return updated_set
        
        # If the list of outliers is empty, then simply return the set.
        else:
            return updated_set

    # The minimize_dfa method recursively splits the initial sets until they no longer change
    def minimize_dfa(self, dfa, split_set, alphabet, prev=None):
        # Save the set before it is modified
        current_set = split_set

        # Go over every set in the list
        for split_states in split_set:
            # If there is only 1 element, don't check
            if len(split_states) < 2:
                pass

            else:
                # Otherwise first check the transitions, identify outliers within the current set and split them
                transition_list = self.check_transitions(dfa, split_set, split_states, alphabet)
                minimized_set = self.identify_outliers(transition_list, split_set)
                # Update the input set
                split_set = minimized_set

        # Check the saved set against the updated set, if it changed run the function again
        if current_set != minimized_set:
            return self.minimize_dfa(dfa, minimized_set, alphabet, current_set)
        
        # Once there are no more changes, return the minimized dfa set
        else:
            return minimized_set

    # The get_transition_table method matches the new combined states and get the union of their transitions
    def get_transition_table(self, input_dfa, minimized_dfa_sets, symbols):
        minimized_dfa = {}

        # Go over every set of combined states
        for state_set in minimized_dfa_sets:

            # If ther eis only a single state, then simply copy over the transitions
            if len(state_set) < 2:
                minimized_dfa[state_set[0]] = input_dfa[state_set[0]]

            else:
                # Otherwise create a new list to keep track of symbols of the states' transitions
                new_transitions = [[] for i in range(len(symbols))]

                # For each state get their transitions
                for state in state_set:
                    state_transitions = input_dfa[state]
                    
                    # Add them to the new list depending on their transition symbol
                    for index, transition in enumerate(state_transitions):
                        if len(new_transitions[index]) == 0 or transition not in new_transitions[index]:
                            new_transitions[index].append(transition)

                # Separate transitions that are singular
                for index in range(len(new_transitions)):
                    if len(new_transitions[index]) == 1:
                        new_transitions[index] = new_transitions[index][0]
                        
                # Add this unon to the newly made dictionary
                minimized_dfa[str(state_set)] = new_transitions

        return minimized_dfa
    
            
if __name__ == "__main__":
    # Change alphabet when more characters
    symbols = [0, 1]

    # Initial state indicated by -
    # Final state indicated by *
    # Input format follows alphabet
    input_dfa = {"-q0": ["q1", "q2"],
                 "q1": ["q1", "q3"],
                 "q2": ["q1", "q2"],
                 "q3": ["q1", "*q4"],
                 "*q4": ["q1", "q2"]
                 }

    input_dfa2 = {
                    "-A": ["B", "F"],
                    "B": ["G", "*C"],
                    "*C": ["*C", "-A"],
                    "D": ["*C", "G"],
                    "E": ["H", "F"],
                    "F": ["*C", "G"],
                    "G": ["G", "E"],
                    "H": ["G", "*C"]
                }

    DFA_Minimizer = DFA_Equivalence_Minimizer(input_dfa2, symbols)    
    print("The minimized DFA set:")
    print(DFA_Minimizer.get_minimized_set(), "\n")

    print("New Transition Table:")
    minimized_dfa = DFA_Minimizer.get_minimized_dfa()
    for state, transitions in minimized_dfa.items():
            print(state, ":", transitions)

