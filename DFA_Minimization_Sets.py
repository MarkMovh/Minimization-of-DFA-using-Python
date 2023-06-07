import copy

# Notes:
# 1. Finish adding comments
# 2. Fix notation in "identify_outliers" function
# 3. Add function to see minimized DFA transition table

class InvalidDeterministicFiniteAutomata(Exception):
    "Invalid Finite Automata Input Detected"
    pass

# Function ascertains that the input fa is a valid deterministic automata
def check_dfa(fa, symbols):
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

        return initial_split(fa, symbols)

    except InvalidDeterministicFiniteAutomata:
        print('''
Error: The following automata you have input is not deterministic.

An Automata can only be minimized if it is in deterministic form.
Convert the input automata into DFA first.
            ''')

# Get all the states and separate them into two sets:
# - One containing the final states
# - One containing all other states
def initial_split(dfa, alphabet):
    initial_set = [[], []]

    # The split is done by checking the initial of the state for the symbol *
    for state in dfa.keys():
        if state[0] == "*":
            initial_set[0].append(state)
        else:
            initial_set[1].append(state)

    return minimize_dfa(dfa, initial_set, alphabet)  

# The transitions of each state are checked, and summarised in a dictionary
def check_transitions(dfa, full_set, inner_set, alphabet):
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

# The following function determines which states should be moved from their list.
def identify_outliers(dfa, transitions, previous_set, alphabet):
    count_dict = {}

    # Iterate over every state and its list locations, saving them to the dictionary
    for value in transitions.values():
        count_dict[tuple(value)] = count_dict.get(tuple(value), 0) + 1

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

def minimize_dfa(dfa, split_set, alphabet, prev=None):
    current_set = split_set

    for split_states in split_set:
        if len(split_states) < 2:
            pass

        else:
            transition_list = check_transitions(dfa, split_set, split_states, alphabet)
            minimized_set = identify_outliers(dfa, transition_list, split_set, alphabet)
            split_set = minimized_set

    if current_set != minimized_set:
        return minimize_dfa(dfa, minimized_set, alphabet, current_set)
    else:
        return minimized_set
            
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

    # Print out the new minimized DFA sets
    print("The minimized DFA:")
    print(check_dfa(input_dfa2, symbols), "\n")

    # Print out the new transition table of Minimized DFA
    print("New Transition Table:")

