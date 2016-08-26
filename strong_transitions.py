def read_transitions(filename):
    trans_strength = []  # An array to hold the gaps and associated transition strengths
    f = open(filename,'r')
    for file_length, l in enumerate(f.readlines()):
        line = l.strip()
        line = line.split()
        if len(line) == 26:
            trans_strength.append([float(line[17]),float(line[20])])
    return trans_strength

def strong_transitions(data, threshold):
    ''' Return the transitions with a dipole matrix element grater than threshold'''
    st = []
    at = []
    for element in data:
        if element[1] > threshold:
            st.append(element[0])
        at.append(element[0])
    return st,at
    
ts = read_transitions('Transitions.dat')
strong_transitions, all_transitions = strong_transitions(ts,0.1)

print "Smallest strong transition:", min(strong_transitions)



