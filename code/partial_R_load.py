import sys
import os
import time

if __name__=='__main__':

    time0 = time.time()

    DESC_FILE = 'n1024/desc/seed1774555967/desc.total.rels.CHKPT.714'

    uv_primes = [
        2,
        2,
        983,
        1511,
        13007,
        101879,
        359449,
        2392163,
        3339857,
        480217663,
        781069339,
        2154642821,
        5457314908915093633,
        69493808444443143151,
        2666515570784707379587,
        22987445119726781943839,
        401483176811778201414487,
        61,
        1009,
        4729,
        240719,
        318403,
        5717603,
        206768819,
        454612757,
        860861137,
        3185039960201617373,
        31905592327134363217,
        449111919860754455821,
        570917114759400722614787,
        628050395217726952177507
    ]

    desc_primes = set()
    # technically, just the norm. it doesn't matter too much
    # if we load some extras in there.

    for jj in uv_primes:
        if jj < 2**36:
            desc_primes.add(jj)

    with open(DESC_FILE, 'r') as descfile:
        for line in descfile:
            if line.startswith('#'):
                continue

            if 'Taken' in line:
                relation_str = line.strip().split(" ")[1]
                rat_str = relation_str.split(":")[1]
                alg_str = relation_str.split(":")[2]

                for rat_fac in rat_str.split(","):
                    rat_prime = int(rat_fac, 16)
                    if rat_prime < 2**36:
                        desc_primes.add(rat_prime)

                for alg_fac in alg_str.split(","):
                    alg_prime = int(alg_fac, 16)
                    if alg_prime < 2**35:
                        desc_primes.add(alg_prime)

    print(f"Encountered {len(desc_primes)} primes in the descent output.")

    RENUMBER_INFO_FILE = 'n1024/renumber.explained'
    PARTIAL_RENUMBER_FILE = 'n1024/partial.renumber.map'

    # We have a few constants to save:
    #   number_of_rational_queries, number_of_fb_valuations
    # Plus mappings between:
    #   rational_prime_to_prime_index
    #   renumber_to_column (-1 for rational primes)
    #   renumber_to_rational_prime_index
    #   rational_prime_index_to_prime
    #   column_to_renumber
    #   column_to_sage_ideal
    #
    # The J ideal is always column 0 = renumber 0.
    # 2 is always the rational prime at index 0.
    #
    # This file will contain the following info:
    # <renumber index> <side-restricted index> <the entire line from renumber.explained>

    writing_file = open(PARTIAL_RENUMBER_FILE, 'w')
    current_renumber_index = 0
    current_column = 0
    current_rational_prime_index = 0

    num_alg_primes_upto_31 = -1

    with open(RENUMBER_INFO_FILE, 'r') as reading_file:
        for line in reading_file:
            if line.startswith('#'):
                continue

            line = line.strip()
            #renum_index_str = line.split()[0]
            #assert renum_index_str[0:4] == 'i=0x'
            #renum_index = int(renum_index_str[4:], 16)
            #assert renum_index == current_renumber_index

            parser, *data = line.split()

            if parser == 'J':
                # definitely include J
                writing_file.write(f'{current_renumber_index} {current_column} {line}\n')
                current_column += 1

            elif parser == 'rat':
                side, p_str = data
                assert int(side) == 0
                this_prime = int(p_str)

                if this_prime in desc_primes:
                    writing_file.write(
                        f'{current_renumber_index} {current_rational_prime_index} {line}\n'
                    )

                current_rational_prime_index += 1

            else:
                assert parser in ['easy', 'proj', 'generic']
                side_str = line.split()[1]
                assert side_str == '1'
                p_str = line.split()[2]
                this_prime = int(p_str)

                if this_prime in desc_primes:
                    writing_file.write(f'{current_renumber_index} {current_column} {line}\n')

                current_column += 1

                if num_alg_primes_upto_31 == -1 and this_prime > 2**31:
                    num_alg_primes_upto_31 = current_column

            current_renumber_index += 1

    writing_file.close()

    total_num_rational_primes = current_rational_prime_index
    total_num_alg_primes = current_column

    time1 = time.time()
    print(f"Took time: {time1-time0}")

    print(f"num_alg_primes_upto_31: {num_alg_primes_upto_31}")
    print(f"total_num_rational_primes: {total_num_rational_primes}")
    print(f"total_num_alg_primes: {total_num_alg_primes}")
