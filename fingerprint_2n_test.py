# Ben Williams '25
# benjamin.r.williams.25@dartmouth.edu or roaf676@gmail.com
# Randomized Algorithms, May 23rd, 2023
# Testing the "unique fingerprint experiment" with the largest highly composite numbers known to man

from fingerprinting_2n_hash_experiment import input_ab_experiment_deterministic
from fingerprinting_2n_hash_experiment import input_ab_experiment_rand
from fingerprinting_2n_hash_experiment import find_succeeds_2n

import math
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick


start_line = 779174

# For storing the output of the worst HCN from the test


# Reads the primes from the primes1.txt (so I don't have to calculate them)
def get_primes_from_file(file_loc, max_count):
    file = open(file_loc, 'r')
    primes_list = []
    prime_count = 0
    line_count = 0
    for line in file:
        if line_count % 2 != 1:
            line = line.strip()
            line = line.strip("\n")
            nums = line.split()

            for num in nums:
                primes_list.append(int(num))
                prime_count += 1
        if prime_count > max_count:
            return primes_list

        line_count += 1

    return primes_list


primes = get_primes_from_file("./files/primes1.txt", 100000)


def test_highly_composites_from_file(file_input, file_output, max_line=None, output_all=False, print_info=False, random_test=False, start_line=None, plot_failure=False):
    # Open the files (we assume they are not None)
    hcn_file = open(file_input, 'r')
    output_file = open(file_output, "w")

    # Init all the variables
    highly_comps = []
    line_count = 1
    fail_rates = []
    indices = []
    log2_list = []
    max_fail_rate = 0
    worst_hcn = (0, 0, [])

    # Read lines ahead of time and do nothing until we are at the start line
    if start_line:
        while line_count < start_line:
            hcn_file.readline()
            line_count += 1
    else:
        # Otherwise just ignore the first line since its 0
        hcn_file.readline()

    # Loop through every HCN (each line is an HCN)
    for line in hcn_file:
        # Stop if we are at the max-wanted HCN
        if max_line:
            if line_count > max_line:
                break

        # Init hcn
        hcn = 1

        # Handle the line
        line = line.strip()
        line = line.strip("\n")
        num_pairs = line.split()

        # The first entry is the number of unique primes
        num_primes = num_pairs.pop(0)
        relevant_primes = []

        # Make a smaller list (maybe not necessary)
        for i in range(int(num_primes)):
            relevant_primes.append(primes[i])

        # Keep track of which prime we are multiplying
        prime_index = 0
        prev_prime_index = 0
        prime_list = [] # in the format of (prime, exponent)

        # Loop through each set of exponent ^ number of primes it is repeated for
        for pair in num_pairs:
            # first entry is exponent, second is number of times we repeat said exponent
            array = pair.split("^")

            # If we are repeating the exponent for more than one prime...
            if len(array) == 2:
                while prime_index < int(array[1]) + prev_prime_index:
                    hcn *= relevant_primes[prime_index] ** int(array[0])
                    prime_list.append((relevant_primes[prime_index], int(array[0])))
                    prime_index += 1
            else:
                # If we are only using the exponent for this single prime
                hcn *= relevant_primes[prime_index] ** int(array[0])
                prime_list.append((relevant_primes[prime_index], int(array[0])))
                prime_index += 1

            prev_prime_index = prime_index

        highly_comps.append(hcn)

        # Calculate the log_2 of the hcn, so that we know what n to use for the experiment
        log2 = math.log(hcn, 2)
        log2_list.append(log2)

        # Run the experiment, add one to the log2 to correspond with the binary
        # Either run the test with every possible factor 1->2n, or with random r's
        if not random_test:
            fail_rate = input_ab_experiment_deterministic((0, hcn), math.ceil(log2) + 1)
        else:
            fail_rate = input_ab_experiment_rand((0, hcn), math.ceil(log2) + 1, 1000)

        indices.append(line_count)
        fail_rates.append(fail_rate)

        # Keep track of which fail rate is the worst, and its respective hcn
        if fail_rate > max_fail_rate:
            max_fail_rate = fail_rate
            worst_hcn = (max_fail_rate, hcn, prime_list)

        # If we want to print info for every hcn
        if print_info:
            if hcn != 0:
                print(f'log2 : {log2}, hcn: {hcn}')
                print(f'Fail rate for hcn # {line_count} = {fail_rate}')

        # If we want to output only the worst HCN per 1000 (obviously much faster)
        if not output_all:
            if line_count % 1000 == 0:
                output_file.write(f'By HCN # {line_count} with log2 = {log2}, Worst Fail Rate = {worst_hcn[0]}, HCN = {worst_hcn[1]} \n')
                print(f'Line number: {line_count}')
        else:
            # or output every HCN's fail rate to the file
            output_file.write(str((fail_rate, hcn)) + "\n")

        line_count += 1

    log2_worst = math.ceil(math.log(worst_hcn[1], 2))
    print(f'\nWorst hcn: {worst_hcn[1]} \nLog2: {log2_worst} \nFail Rate: {worst_hcn[0]} \nPrimes and Exponents:{worst_hcn[2]}')

    # Plot the probability of failure for all the HCNs we have seen
    if plot_failure:
        percent_fail_rates = []
        for fr in fail_rates:
            percent_fail_rates.append(fr * 100)
        plot = plt.subplot2grid((2, 2), (0, 0), rowspan=2, colspan=2)
        plot.plot(log2_list, percent_fail_rates)

        # Set the title of the plot accordingly
        if start_line:
            if max_line:
                plot.set_title(f'Probability of failure from HCN #{start_line} to HCN #{max_line}')
            else:
                plot.set_title(f'Probability of failure from HCN #{start_line} to HCN #779674')
        elif max_line:
            plot.set_title(f'Probability of failure from the first HCN to HCN #{max_line}')
        else:
            plot.set_title(f'Probability of failure for HCNs as n grows')

        plot.set_xlabel("Log_2(HCN)")
        plot.yaxis.set_major_formatter(mtick.PercentFormatter())
        plot.set_ylabel("Percentage of numbers <= 2n as factors")
        plt.show()

    # Only doable for smaller HCNs, since they just get too big to plot and matplotlib crashes
    # plt.plot(highly_comps, fail_rates)
    # plt.show()

    return worst_hcn[1], log2_worst, worst_hcn[2]

# Test the first 10,000 HCNs
# test_highly_composites_from_file("./files/HCN.txt", "./files/null.txt", max_line=10000, plot_failure=True)
# test_highly_composites_from_file("./files/HCN.txt", "./files/null.txt", start_line=start_line, plot_failure=True)


# Use the information to optimize the HCN
# I changed the start line and max line to closely straddle the worst HCN in the last 500 for my own sanity
triple = test_highly_composites_from_file("./files/HCN.txt", "./outputs/null.txt", start_line=779374, max_line=779474)


# Recalculates a number given its list of primes and their exponents
# Needs to be done since we crash when dividing along the way
def recalc_candidate(primes_and_exponents):
    result = 1
    for pe in primes_and_exponents:
        result *= (pe[0] ** pe[1])
    return result


# Takes any number as the parameter and gets its prime factorization
def get_prime_factors(num, print_info=False):
    copy = num
    i = 2
    index_map = dict()
    prime_factors = []
    if print_info:
        print(f'Getting prime factors for {num}')

    while i * i < copy:
        while copy % i == 0:
            if index_map.get(i) is None:
                prime_factors.append((i, 1))
                index_map[i] = len(prime_factors) - 1
            else:
                prime_factors[index_map[i]] = (i, prime_factors[index_map[i]][1] + 1)
            copy /= i
        i += 1

    if len(prime_factors) == 0:
        prime_factors.append((i, 1))
        return prime_factors, None

    return prime_factors, index_map


# Counts the number of primes from a start prime and a finish number
# (here, the last prime in the prime factorization, and 2*n as the finish)
def missing_primes_count(start, finish):
    count = 0
    # Define this to be the total number of primes and multiples of primes that we are missing below 2n
    total_wrong_count = 0

    # Find the start prime in the list
    index = 0
    while primes[index] != start:
        index += 1

    # to account for the "start" prime
    index += 1
    while primes[index] < finish:
        total_wrong_count += 1
        count += 1

        mult = 2
        while primes[index] * mult < finish:
            total_wrong_count += 1
            mult += 1
        index += 1

    return count, total_wrong_count


# Takes a HCN as input, prints out the info about the best
def optimize_hcn_for_2n_attempt(hcn, log2, primes_and_exponents):
    log2 = math.ceil(log2)
    old_fail_rate = input_ab_experiment_deterministic((0, hcn), log2)
    succeeding_nums = find_succeeds_2n((0, hcn), log2)
    candidate_pe = primes_and_exponents

    remove_count = 0
    # Get rid of all excess exponents
    for i in range(len(candidate_pe)):
        while candidate_pe[i][0] ** candidate_pe[i][1] > 2*log2:
            candidate_pe[i] = (candidate_pe[i][0], candidate_pe[i][1] - 1)
            remove_count += 1
    print(f'Removed {remove_count} prime-exponents')

    candidate = recalc_candidate(candidate_pe)
    #print(f'Old log2: {log2}, New log2: {math.log(candidate, 2)}')
    diff_count = 0

    # Takes in the first number that succeeded (also the smallest)
    remove = succeeding_nums.pop(0)
    remove_prime_factors, index_map = get_prime_factors(remove)

    while candidate * remove < 2**log2:
        # Find the differences in prime factorization, and just multiply those, rather than the full "remove":
        prime_differences = []
        difference = False
        for i in range(len(remove_prime_factors)):
            for j in range(i, len(candidate_pe)):
                # If we are at the same prime factor
                if remove_prime_factors[i][0] == candidate_pe[j][0]:
                    if remove_prime_factors[i][1] > candidate_pe[j][1]:
                        prime_differences.append((remove_prime_factors[i][0], remove_prime_factors[i][1] - candidate_pe[j][1]))
                        difference = True

        # The removed number is a prime
        if not difference:
            candidate_pe.append((remove, 1))
            candidate *= remove
            diff_count += 1
        else:  # Otherwise add the necessary different prime factors between
            for diff in prime_differences:
                for i in range(len(candidate_pe)):
                    if candidate_pe[i][0] == diff[0]:
                        candidate_pe[i][1] += diff[1]
                candidate *= (diff[0] ** diff[1])
                diff_count += diff[1]

        # More than one number could have been removed by adding these prime
        succeeding_nums = find_succeeds_2n((0, candidate), log2)
        # print(f'Current Fail Rate: {1 - (len(succeeding_nums) / (2*log2))}')

        if len(succeeding_nums) == 0: # If by some miracle we succeed...
            print(f'Success with {candidate}')
        else:
            remove = succeeding_nums.pop(0)
            remove_prime_factors, index_map = get_prime_factors(remove)

    new_fail_rate = 1 - (len(succeeding_nums) / (2*log2))
    print(f'Total prime differences = {diff_count}')
    print(f'Best optimized fail rate: {1 - (len(succeeding_nums) / (2*log2))}\nNumber: {candidate}\nImproved(?) by: {new_fail_rate - old_fail_rate}\n New prime factorization:\n{candidate_pe}')

    # To figure out how many primes we are missing before 2*n
    final_prime = candidate_pe[len(candidate_pe) - 1][0]
    missing_count = missing_primes_count(final_prime, 2*log2)
    print(f'Missing {missing_count[0]} primes before 2n = {2*log2}')
    print(f'Missing {missing_count[1]} primes and prime-multiples before 2n = {2 * log2}')
    return new_fail_rate, candidate


# Only run the final algorithm if we can
if triple:
    output = open("../outputs/optimized_hcn_final_last_500", 'w')
    info = optimize_hcn_for_2n_attempt(triple[0], triple[1], triple[2])
    output.write("Optimized HCN: " + str(info[1]) + "\nFail Rate:" + str(info[0]) + "\nLog2: " + str(math.ceil(math.log(info[1], 2))))







