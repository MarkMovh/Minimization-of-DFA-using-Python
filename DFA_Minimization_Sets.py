import copy

# 1. Fix the final transition table creations to include transition combinations of more than 2
# 2. Make each class inherit methods that are similar / identical
# y. Output a node graph for the minimized DFA
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

    # Get methods for the initialized information
    def get_minimized_set(self):
        print(self.minimized_set)
    
    # Get method to show the minimized dfa
    def get_minimized_dfa(self):
        print("New Transition Table:")
        for state, transitions in self.minimized_dfa.items():
                print(state, ":", transitions)


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

            # If there is only a single state, then simply copy over the transitions
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
    
# DFA minimizer method based on Myhill-Nerode method
class DFA_Myhill_Nerode_Minimizer(DFA_Scanner):
    # Initialize the class by getting the minimized set and dfa
    def __init__(self, dfa, alphabet):
        self.dfa = dfa
        self.alphabet = alphabet

        self.check_dfa()

        initial_table = self.initial_min_table(self.dfa)
        self.minimized_table = self.fill_table(self.dfa, initial_table)
        self.minimized_dfa = self.get_transition_table(self.minimized_table, self.dfa, self.alphabet)

    # Get method to show the final minimized table
    def get_minimized_table(self):

        for index, combination in enumerate(self.minimized_table):
            print(list(self.dfa.keys())[index], combination)

    # Get method to show the minimized dfa
    def get_minized_dfa(self):
        print("New Transition Table:")
        for state, transitions in self.minimized_dfa.items():
                print(state, ":", transitions)

    # Using inheritance, a check method can be created to ensure the FA is DFA
    def check_dfa(self):
        return DFA_Scanner.check_dfa(self, self.dfa, self.alphabet)

    # Create an initial table comparing state transitions
    def initial_min_table(self, dfa):
        init_table = []

        # Generate an empty matrix
        for state_num in range(len(dfa.keys())):
            init_table.append([0 for transitions in range(len(dfa.keys()))])

        # Mark off pairs where one is final and other non-final
        for row, state_1 in enumerate(dfa.keys()):
            for column, state_2 in enumerate(dfa.keys()):
                if ("*" in state_1 and "*" not in state_2) or ("*" not in state_1 and "*" in state_2):
                    init_table[row][column] = 1
                else:
                    pass
                
        return init_table

    # Begin filling in the table based off the initial table
    def fill_table(self, dfa, table):
        input_table = copy.deepcopy(table)

        # Check every row and column
        for row, state_1 in enumerate(dfa.keys()):
            for column, state_2 in enumerate(dfa.keys()):
                # Pass if 1, as it is already filled
                if (table[row][column] == 1):
                    continue
                # Pass if each state is a final state
                if "*" in state_1 and "*" in state_2:
                    pass
                # If the states are the same, mark as "x"
                if state_1 == state_2:
                    table[row][column] = "x"
                # If none are true, check the combination transitions
                else:
                    connection = self.check_combination(state_1, state_2, dfa, table)

                    # If transitions are unique, mark as 1
                    if connection:
                        table[row][column] = 1
                    else:
                        pass

        # Recursively fill in table, base case being that there is no table change from the input table.
        if input_table != table:
            return self.fill_table(dfa, table)
        else:
            return table

    # The following method checks if transitions are marked on the current table
    def check_combination(self, state1, state2, dfa, comb_table):
        state1_transitions = dfa[state1]
        state2_transitions = dfa[state2]

        combination = False

        # check each transition combination in the table
        for symbol_num in range(len(state1_transitions)):
            # Get the location based on the symbol transitions
            row = list(dfa.keys()).index(state1_transitions[symbol_num])
            column = list(dfa.keys()).index(state2_transitions[symbol_num])

            # If it is marked as 1, then return true for "to be marked"
            if comb_table[row][column] == 1:
                combination = True
                break

            else:
                pass
        
        return combination
    
    # Method that creates the final transition table
    def get_transition_table(self, minimized_table, input_dfa, symbols):
        minimized_set = self.create_minimized_dfa_set(minimized_table, input_dfa)
        minimized_dfa = {}

        # Go over every set of combined states
        for state_set in minimized_set:

            # If there is only a single state, then simply copy over the transitions
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

    # Method that constructs the minimized dfa set
    def create_minimized_dfa_set(self, min_table, dfa):
        minimized_dfa_set = []

        # Iterate over each row and column, save index of each to avoid duplicate values
        for r_index, row in enumerate(min_table):
            for c_index, col_value in enumerate(row):
                row_state = list(dfa.keys())[r_index]
                col_state = list(dfa.keys())[c_index]
                # Ignore the top right half of the matrix
                if c_index > r_index:
                    break
                # Ignore x and 1 as they are not relevant
                elif col_value == "x" or col_value == 1:
                    pass
                # 0 highlights an unmarked pair, ergo combine them into a single transition
                else:
                    minimized_dfa_set.append([col_state, row_state])

                    # Additionally, check if the matching state is already in another state
                    for transition in minimized_dfa_set:
                        # If so, and the current checking state is not in there, add it and remove the just added combination
                        if col_state in transition and row_state not in transition:
                            transition.append(row_state)
                            minimized_dfa_set.pop()
                        else:
                            pass

                    break

        # Retrieve a flattened set of the matrix
        flattened_set = [state for combination in minimized_dfa_set for state in combination]

        # Check each of the original states to see if they have been added.
        for state in list(dfa.keys()):
            if state not in flattened_set:
                minimized_dfa_set.append([state])

        return minimized_dfa_set



if __name__ == "__main__":
    # Change alphabet when more characters
    symbols = [0, 1]

    # Initial state indicated by -
    # Final state indicated by *
    # Input format follows alphabet
    input_dfa = {
                "-q0": ["q1", "q2"],
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
    
    input_dfa3 = {
                "-q0": ["q3", "*q1"],
                 "*q1": ["*q2", "q5"],
                 "*q2": ["*q2", "q5"],
                 "q3": ["-q0", "*q4"],
                 "*q4": ["*q2", "q5"],
                 "q5": ["q5", "q5"]
    }

    # DFA_Minimizer = DFA_Equivalence_Minimizer(input_dfa3, symbols)    
    # DFA_Minimizer.get_minimized_set()
    # DFA_Minimizer.get_minimized_dfa()
    
    DFA_Minimizer = DFA_Myhill_Nerode_Minimizer(input_dfa3, symbols)
    DFA_Minimizer.get_minimized_table()
    DFA_Minimizer.get_minized_dfa()
