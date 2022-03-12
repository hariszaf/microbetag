#!/usr/bin/env python

import json
import glob
import sys

def flatten(list_of_lists):
    if len(list_of_lists) == 0:
        return list_of_lists
    if isinstance(list_of_lists[0], list):
        return flatten(list_of_lists[0]) + flatten(list_of_lists[1:])
    return list_of_lists[:1] + flatten(list_of_lists[1:])

# def get_modules_of_genome_under_study(genome):
def count_nested_lists(l):
   count = 0 
   for e in l: 
      if isinstance(e, list):
         count = count + 1 + count_nested_lists(e)
   return count

# Remove optional parts if any
def check_for_minus(alternative):

   complex_check = False

   if isinstance(alternative, list):
      tmp = []
      for i in alternative:    

         if "-" in i:
            tmp.append(i.split("-")[0])

         else:
            tmp.append(i)
      alternative = tmp

   else: 
      if "-" in alternative: 
         alternative = alternative.split("-")[0]

   return alternative

# Get missing parts from a species modules
def extract_missing_parts(ncbi_id_beneficiary):
  
   f      = open("../../ref-dbs/module_definition_map.json", "r")
   mo_map = json.load(f)

   for file in glob.glob("../../ref-dbs/kegg_genomes/" + str(ncbi_id_beneficiary) + "/*.json"):
      with open(file, 'r') as f:
         genome_mos = json.load(f)

   # Species in-need of by-products
   count_of_species_modules = 0
   parts_missing_from_modules = {}
   beneficary_kos_per_module = {}

   # Parse genome's KEGG modules and terms as found in the KEGG ORGANISMS db 
   for module, kos_on_its_own in genome_mos.items(): 

      list_of_kos_present                           = [k[3:] for k in kos_on_its_own.values()]
      beneficary_kos_per_module[module]             = list_of_kos_present
      count_of_species_modules                     += 1
      parts_missing_from_modules[module] = {}

      if module != "md:M00051":
         continue

      for index, step in enumerate(mo_map[module]['steps']):

         step_complet = False
         nested = count_nested_lists(step)

         """
         We will have a list (cases) for the various alternatives that will
         come up as sets of KOs that if found can complete the step
         We will only add cases in the parts_missing_from_modules 
         dictionary if and oly if there is not already a way for a complete step
         """
         parts_missing_from_modules[module][index] = {}
         cases = []

         """
         In this case the children of the step, are the alternatives 
         and for each alternative, all the KOs are needed for the step
         to be complete
         """
         if nested > 0: 

            step_alternatives_map = {}

            for option, alternative in enumerate(step): 

               if isinstance(alternative, str):
                  alternative = [alternative]

               # Check for minus or/and pluses in step
               # print("-------- alternative: ", alternative)
               alternative_to_check = check_for_minus(alternative)
               # print("========= alternative to check: ", alternative_to_check)
               checking = alternative_to_check[:]
               # print("********** checking: ", checking)

               # Investigate the alternatives
               for inner_index, ko_of_alternative in enumerate(alternative_to_check): 

                  if "+" in ko_of_alternative: 

                     complex   = ko_of_alternative.split("+")
                     compounds = complex[:]

                     for term in complex: 

                        if term in list_of_kos_present:
                           compounds.remove(term)
                     
                     if len(compounds) == 0:
                        checking.remove(ko_of_alternative)

                     else:
                        checking.remove(ko_of_alternative)
                        checking.append(compounds)

                  if " " in ko_of_alternative: 

                     multistep = ko_of_alternative.split(" ")
                     compounds = multistep[:]

                     for term in multistep: 

                        if term in list_of_kos_present:
                           compounds.remove(term)
                     
                     if len(compounds) == 0:
                        checking.remove(ko_of_alternative)

                     else:
                        checking.remove(ko_of_alternative)
                        checking.append(compounds)

                  else:

                     """
                     check if something's missing...
                     """

                     if ko_of_alternative in list_of_kos_present:
                        checking.remove(ko_of_alternative)

               if len(checking) == 0: 
                  step_complet = True

               else: 
                  step_alternatives_map[option] = checking

            if step_complet == False: 
               parts_missing_from_modules[module][index] = step_alternatives_map

         else:
            """
            In that case the alternatives, if any, are at the top level 
            and only one of them is needed for the step to be complete
            """

            step_without_signs = check_for_minus(step)
            step = step_without_signs[:]

            if isinstance(step, str):
               """
               The step can be as a string
               """
               if "+" in step: 

                  complex       = step.split("+")
                  complex_check = complex[:]
                  # print(complex) is missi

                  for ko in complex:
                     if ko in list_of_kos_present:
                        complex_check.remove(ko)

                  if len(complex_check) == 0:
                     # print("Complex complete on the species")
                     continue

                  else:
                     parts_missing_from_modules[module][index] = [complex_check]

            else:

               complexes = False
               tmps      = {}

               for case, alternative in enumerate(step):

                  if "+" in alternative: 
                     """
                     This is the case where a complex is on the step
                     """
                     complexes     = True
                     complex       = alternative.split("+")
                     complex_check = complex[:]

                     for ko in complex:
                        if ko in list_of_kos_present:
                           complex_check.remove(ko)

                     if len(complex_check) == 0:
                        tmps[case] = []

                     else:
                        # parts_missing_from_modules[module][index].append(complex_check)
                        tmps[case] = complex_check

                  else:
                     """
                     This is the case where there is no complex on the step
                     """

                     if alternative in list_of_kos_present:
                        step_complet = True
                        parts_missing_from_modules[module][index] = []
                        break

               if step_complet == False: 

                  if complexes == False:
                     for case, i in enumerate(step):
                        parts_missing_from_modules[module][index][case] = [i]
                  else:
                     if len(tmps) == 0:
                        for case, i in enumerate(step):
                           parts_missing_from_modules[module][index][case] = [i]
                     else:
                        for case, i in enumerate(step):
                           if "+" in i:
                              parts_missing_from_modules[module][index][case] = tmps[case]
                           else:
                              parts_missing_from_modules[module][index][case] = [i]

   missing_steps = {}
   for module, missings in parts_missing_from_modules.items():       
      if missings: 
         missing_steps[module] = {}
         list_of_kos_present   = beneficary_kos_per_module[module]
         for step, cases in missings.items(): 
            missing_steps[module][step] = {}
            if isinstance(cases, list):
               flat_case = flatten(cases)
               flat_case = [x for x in flat_case if x not in list_of_kos_present]
               missing_steps[module][step]['0'] = flat_case
            elif isinstance(cases, dict):
               for number, case in cases.items():
                  flat_case = flatten(case)
                  flat_case = [x for x in flat_case if x not in list_of_kos_present]
                  missing_steps[module][step][number] = flat_case

   return missing_steps

# Get complementary cases if any from its partner
def check_for_complements(ncbi_id_donor, missing_parts, beneficiary_ncbi_id):

   for file in glob.glob("../../ref-dbs/kegg_genomes/" + str(ncbi_id_donor) + "/*.json"):
      with open(file, 'r') as f:
         donors_genome_mos = json.load(f)

   metabolic_interactions_per_module = {}

   for module, lacking_steps in missing_parts.items():

      print("MODULE ", module)

      if module not in donors_genome_mos.keys():
         continue

      else:

         doners_module = list(donors_genome_mos[module].values())
         doners_module = [x[3:] for x in doners_module]

         print("\n > DONERS MODULE: ", doners_module, "MD: ", module)
         print("\n > MISSING STEPS: ", lacking_steps)

         for beneficiarys_step_index, kos_missing in lacking_steps.items():

            """
            Remember! kos_missing might have 1 or more alternatives
            BUT all the terms of an alternative are required to 
            have a step completed // for now it might be a list or a dict
            I need to fix that... make all dicts
            """

            for missing_ko_alternative in kos_missing.values(): 
               if len(missing_ko_alternative) > 0:
                  check = missing_ko_alternative[:]
                  for potential_interaction in doners_module:
                     if potential_interaction in check:
                        print(" HERE IS A POSSIBLE INTERACTION!", potential_interaction)
                        check.remove(potential_interaction)

                  if len(check) == 0:
                     print("USING THIS WAY, THE MODULE IS NOW COMPLETE!")


   # if len(list_of_metabolic_interactions) > 0: 
   #    print("Species with NCBI Taxonomy Id: ", beneficiary_ncbi_id, "gets KOs", list_of_metabolic_interactions, " from species with NCBI Taxonomy Id: ", doner_ncbi_id)


def main(beneficiary_ncbi_id, doner_ncbi_id):

   missing_parts     = extract_missing_parts(beneficiary_ncbi_id)
   complementarities = check_for_complements(doner_ncbi_id, missing_parts, beneficiary_ncbi_id)
   return complementarities


if __name__=="__main__":

   """
   35128   : tps
   1891787 : rpon
   562     : eco
   """
   main(1891787, 562)
   #print(to_find)

