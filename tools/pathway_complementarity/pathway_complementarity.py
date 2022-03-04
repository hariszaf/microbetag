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

def extract(lists, d):

   output = []

   if d == 1:
      return output.extend(lists)

   for sub_list in lists:
      extract(sub_list, d - 1)

def check_for_signs(alternative):
   if isinstance(alternative, list):
      tmp = []
      for i in alternative:    
         if "-" in i:
            tmp.append(i.split("-")[0])
         elif "+" in i:
            split = i.split("+")
            for j in split:
               tmp.append(j)
         else:
            tmp.append(i)
      alternative = tmp
   else: 
      if "-" in alternative: 
         alternative = alternative.split("-")[0]
      if "+" in alternative: 
         alternative = alternative.split("+")
   return alternative



def main():
  
   f      = open("../../ref-dbs/module_definition_map.json", "r")
   mo_map = json.load(f)


   ncbi_id = str(435)
   for file in glob.glob("../../ref-dbs/kegg_genomes/" + ncbi_id + "/*.json"):
      with open(file, 'r') as f:
         genome_mos = json.load(f)

   # Species in-need of by-products
   count_of_species_modules = 0
   missing_steps_from_module_under_study = {}



   for module, kos_on_its_own in genome_mos.items(): 

      count_of_species_modules += 1
      missing_steps_from_module_under_study[module] = {}

      if module != "md:M00154": 
         continue

      for index, step in enumerate(mo_map[module]['steps']):

         step_complet = False

         nested = count_nested_lists(step)

         """
         We will have a list (cases) for the various alternatives that will
         come up as sets of KOs that if found can complete the step
         We will only add cases in the missing_steps_from_module_under_study 
         dictionary if and oly if there is not already a way for a complete step
         """

         missing_steps_from_module_under_study[module][index] = []
         cases = []

         # print("\n\n>>>>>> ", step, index)

         if nested > 0: 
            """
            In this case the children of the step, are the alternatives 
            and for each alternative, all the KOs are needed for the step
            to be complete
            """

            for option, alternative in enumerate(step): 

               # Check for minus or/and pluses in step
               alternative = check_for_signs(alternative)

               # Investigate the alternatives
               alternative_to_check = alternative[:]

               for ko in kos_on_its_own.values():

                  ko = ko[3:]

                  if ko in alternative_to_check: 
                     
                     if isinstance(alternative_to_check, list):
                        alternative_to_check.remove(ko) 

                     else: 
                        # KEEP IN MIND TO CHECK IF THERE IS A CASE OF '+' IN THESE TERMS
                        print("???! --> ", alternative_to_check)
                        step_complet = True
                        break

               if len(alternative_to_check) == 0:
                  
                  step_complet = True
                  print("Step is complete by the species on its own")
                  break

               else: 
                  cases.append(alternative_to_check)

            if alternative_to_check == False:
               missing_steps_from_module_under_study[module][index].append(cases)



         else:
            """
            In that case the alternatives, if any, are at the top level 
            and only on of them is needed for the step to be complete
            """

            print(step)
            step_without_signs = check_for_signs(step)
            print(step_without_signs)

            for ko in kos_on_its_own.values(): 

               ko = ko[3:]

               if ko in step: 

                  print("Step is already present from the species itself")
                  step_complet = True
                  break

            if step_complet == False: 

               missing_steps_from_module_under_study[module][index] = step

   missing_steps_from_module_under_study = {k: v for k, v in missing_steps_from_module_under_study.items() if v}
   return missing_steps_from_module_under_study





if __name__=="__main__":
   
   to_find = main()
   print(to_find)







# print("\n~~~~~\n")
# print(genome_mos['md:M00034'])
# print(missing_steps_from_module_under_study['md:M00034'])
# print("\n~~~~~\n")
# print("From its own: ", genome_mos['md:M00854'])
# print("To get from elsewher: ", missing_steps_from_module_under_study['md:M00854'])

