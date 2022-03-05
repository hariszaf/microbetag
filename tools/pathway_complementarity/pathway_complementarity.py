#!/usr/bin/env python

import json
import glob


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

# Main function
def main(ncbi_id):
  
   f      = open("../../ref-dbs/module_definition_map.json", "r")
   mo_map = json.load(f)

   for file in glob.glob("../../ref-dbs/kegg_genomes/" + ncbi_id + "/*.json"):
      with open(file, 'r') as f:
         genome_mos = json.load(f)

   # Species in-need of by-products
   count_of_species_modules = 0
   missing_steps_from_module_under_study = {}


   # Parse genome's KEGG modules and terms as found in the KEGG ORGANISMS db 
   for module, kos_on_its_own in genome_mos.items(): 

      list_of_kos_present                           = [k[3:] for k in kos_on_its_own.values()]
      count_of_species_modules                     += 1
      missing_steps_from_module_under_study[module] = {}


      print("MODULE: ", module, "\nON MY OWN: ", kos_on_its_own, "\n")
      print("The complete module: ", mo_map[module]['steps'], "\n")

      for index, step in enumerate(mo_map[module]['steps']):

         step_complet = False
         nested = count_nested_lists(step)

         """
         We will have a list (cases) for the various alternatives that will
         come up as sets of KOs that if found can complete the step
         We will only add cases in the missing_steps_from_module_under_study 
         dictionary if and oly if there is not already a way for a complete step
         """
         missing_steps_from_module_under_study[module][index] = {}
         cases = []


         """
         In this case the children of the step, are the alternatives 
         and for each alternative, all the KOs are needed for the step
         to be complete
         """
         if nested > 0: 

            print("STEP: ", step)
            print("KOs present: ", list_of_kos_present)

            step_alternatives_map = {}

            for option, alternative in enumerate(step): 

               if isinstance(alternative, str):
                  alternative = [alternative]

               # Check for minus or/and pluses in step
               alternative_to_check = check_for_minus(alternative)
               checking = alternative_to_check[:]

               # Investigate the alternatives
               for inner_index, ko in enumerate(alternative_to_check): 

                  if "+" in ko: 

                     complex   = ko.split("+")
                     compounds = complex[:]

                     for term in complex: 

                        if term in list_of_kos_present:
                           compounds.remove(term)
                     
                     if len(compounds) == 0:
                        checking.remove(ko)

                     else:
                        checking.remove(ko)
                        checking.append(compounds)

                  else:

                     if ko in list_of_kos_present:
                        checking.remove(ko)

               if len(checking) == 0: 

                  # print("\n>>>> ALL THE STEP IS HERE !!\n")
                  step_complet = True

               else: 
                  print("CHECKING TO ADD ", checking)
                  step_alternatives_map[option] = checking

            if step_complet == False: 
               
               missing_steps_from_module_under_study[module][index] = step_alternatives_map

         else:
            """
            In that case the alternatives, if any, are at the top level 
            and only one of them is needed for the step to be complete
            """

            step_without_signs = check_for_minus(step)
            step = step_without_signs[:]
            print("This is step: ", step)

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
                     print("Complex complete on the species")

                  else:
                     missing_steps_from_module_under_study[module][index] = [complex_check]

            else:

               complexes = False
               tmps      = {}

               for case, alternative in enumerate(step):

                  if "+" in alternative: 
                     """
                     This is the case where a complex is on the step
                     """
                     print("Here is a complex:")

                     complexes     = True
                     complex       = alternative.split("+")
                     complex_check = complex[:]

                     for ko in complex:
                        if ko in list_of_kos_present:
                           complex_check.remove(ko)

                     if len(complex_check) == 0:
                        tmps[case] = []

                     else:
                        # missing_steps_from_module_under_study[module][index].append(complex_check)
                        tmps[case] = complex_check
                        print("~~~", tmps, case, alternative)

                  else:
                     """
                     This is the case where there is no complex on the step
                     """

                     if alternative in list_of_kos_present:
                        step_complet = True
                        missing_steps_from_module_under_study[module][index] = []
                        break

               if step_complet == False: 

                  if complexes == False:
                     for case, i in enumerate(step):
                        missing_steps_from_module_under_study[module][index][case] = [i]
                  else:
                     if len(tmps) == 0:
                        for case, i in enumerate(step):
                           missing_steps_from_module_under_study[module][index][case] = [i]
                     else:
                        for case, i in enumerate(step):
                           print(">case: ", case)
                           print("TMPS:", tmps)
                           if "+" in i:
                              print("-------------", tmps[case])
                              missing_steps_from_module_under_study[module][index][case] = tmps[case]
                           else:
                              missing_steps_from_module_under_study[module][index][case] = [i]


            print("\n\n~~~~\n")

   missing_steps_from_module_under_study = {k: v for k, v in missing_steps_from_module_under_study.items() if v}
   return missing_steps_from_module_under_study





if __name__=="__main__":
   ncbi_id = str(435)
   to_find = main(ncbi_id)
   print(to_find)

