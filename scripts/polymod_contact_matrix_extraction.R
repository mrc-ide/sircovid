# Loading In Required Library Packages
library(socialmixr); library(tidyverse)
polymod <- socialmixr::polymod
participants <- polymod$participants %>%
  filter(country == "United Kingdom") %>%
  select(part_id, part_age) %>%
  mutate(part_age = 5 * floor(part_age/5))
UK_ids <- participants$part_id # Identify UK participant IDs
contacts <- polymod$contacts
participants_by_age <- table(participants$part_age) # Get number of participants by age
number_people_df <- data.frame(age_group = as.numeric(names(participants_by_age)),
                               n = as.vector(unname(participants_by_age))) # Create df of age groups and # participants

# Generate Overall Polymod Matrix
# Processing includes calculating mean age for those where age only available as range
# and assigning to 75 all individuals older than 75 (as we have a 75+ category)
overall_contacts <- polymod$contacts %>%
  filter(part_id %in% UK_ids, !(is.na(cnt_age_exact) & 
                                  !is.na(cnt_age_est_min) & 
                                  is.na(cnt_age_est_max))) %>% 
  mutate(age_range_mean = (cnt_age_est_min + cnt_age_est_max)/ 2, 
         age_group = ifelse(is.na(cnt_age_exact), 
                            ifelse(cnt_age_est_min < 75, 5 * floor(age_range_mean/5), 75), 
                            ifelse(cnt_age_exact < 75, 5 * floor(cnt_age_exact/5), 75))) %>% 
  select(part_id, age_group, everything(), -cnt_age_exact, 
         -cnt_age_est_min, -cnt_age_est_max, -age_range_mean,
         -cont_id, -cnt_gender, -frequency_multi, -phys_contact, -duration_multi) %>% 
  left_join(participants, by = "part_id") %>% # Adding in the age of the participant
  group_by(part_age, age_group) %>% # Grouping by participant age and age of person contacted
  summarise(total = sum(cnt_school + cnt_transport + cnt_work + cnt_home + cnt_otherplace + cnt_leisure)) %>% # Summing over all contacts
  left_join(number_people_df, by = c("part_age" = "age_group")) %>% # adding in the number of people surveyed in each age group
  mutate(per_cap = total/n) %>% # dividing total contacts by number surveyed in each age group
  select(part_age, age_group, per_cap) %>%
  spread(age_group, per_cap) # creating the matrix 
ages <- overall_contacts$part_age 
overall_contacts <- as.matrix(overall_contacts[, -1]) # remove column of age group labels
overall_contacts[is.na(overall_contacts)] <- 0 # change the single NA at element [16, 14] which is an NA as 0 contacts between those 2 age groups were reported
row.names(overall_contacts) <- ages
saveRDS(overall_contacts, file = "scripts/overall_contact_matrix.rds")

# Generate Home Matrix
home_contacts <- polymod$contacts %>%
  filter(part_id %in% UK_ids, !(is.na(cnt_age_exact) & 
                                  !is.na(cnt_age_est_min) & 
                                  is.na(cnt_age_est_max))) %>% 
  mutate(age_range_mean = (cnt_age_est_min + cnt_age_est_max)/ 2, 
         age_group = ifelse(is.na(cnt_age_exact), 
                            ifelse(cnt_age_est_min < 75, 5 * floor(age_range_mean/5), 75), 
                            ifelse(cnt_age_exact < 75, 5 * floor(cnt_age_exact/5), 75))) %>% 
  select(part_id, age_group, everything(), -cnt_age_exact, 
         -cnt_age_est_min, -cnt_age_est_max, -age_range_mean,
         -cont_id, -cnt_gender, -frequency_multi, -phys_contact, -duration_multi) %>% 
  left_join(participants, by = "part_id") %>% # Adding in the age of the participant
  group_by(part_age, age_group) %>% # Grouping by participant age and age of person contacted
  summarise(total = sum(cnt_home)) %>% # Summing over all contacts
  left_join(number_people_df, by = c("part_age" = "age_group")) %>% # adding in the number of people surveyed in each age group
  mutate(per_cap = total/n) %>% # dividing total contacts by number surveyed in each age group
  select(part_age, age_group, per_cap) %>%
  spread(age_group, per_cap) # creating the matrix 
ages <- home_contacts$part_age 
home_contacts <- as.matrix(home_contacts[, -1]) # remove column of age group labels
home_contacts[is.na(home_contacts)] <- 0 # change the single NA at element [16, 14] which is an NA as 0 contacts between those 2 age groups were reported
row.names(home_contacts) <- ages
heatmap(home_contacts, Rowv = NA, Colv = NA)
saveRDS(home_contacts, file = "scripts/home_contact_matrix.rds")

# Generate School Matrix
school_contacts <- polymod$contacts %>%
  filter(part_id %in% UK_ids, !(is.na(cnt_age_exact) & 
                                  !is.na(cnt_age_est_min) & 
                                  is.na(cnt_age_est_max))) %>% 
  mutate(age_range_mean = (cnt_age_est_min + cnt_age_est_max)/ 2, 
         age_group = ifelse(is.na(cnt_age_exact), 
                            ifelse(cnt_age_est_min < 75, 5 * floor(age_range_mean/5), 75), 
                            ifelse(cnt_age_exact < 75, 5 * floor(cnt_age_exact/5), 75))) %>% 
  select(part_id, age_group, everything(), -cnt_age_exact, 
         -cnt_age_est_min, -cnt_age_est_max, -age_range_mean,
         -cont_id, -cnt_gender, -frequency_multi, -phys_contact, -duration_multi) %>% 
  left_join(participants, by = "part_id") %>% # Adding in the age of the participant
  group_by(part_age, age_group) %>% # Grouping by participant age and age of person contacted
  summarise(total = sum(cnt_school)) %>% # Summing over all contacts
  left_join(number_people_df, by = c("part_age" = "age_group")) %>% # adding in the number of people surveyed in each age group
  mutate(per_cap = total/n) %>% # dividing total contacts by number surveyed in each age group
  select(part_age, age_group, per_cap) %>%
  spread(age_group, per_cap) # creating the matrix 
ages <- school_contacts$part_age 
school_contacts <- as.matrix(school_contacts[, -1]) # remove column of age group labels
school_contacts[is.na(school_contacts)] <- 0 # change the single NA at element [16, 14] which is an NA as 0 contacts between those 2 age groups were reported
row.names(school_contacts) <- ages
heatmap(school_contacts, Rowv = NA, Colv = NA)
saveRDS(school_contacts, file = "scripts/school_contact_matrix.rds")

# Generate Work Matrix
work_contacts <- polymod$contacts %>%
  filter(part_id %in% UK_ids, !(is.na(cnt_age_exact) & 
                                  !is.na(cnt_age_est_min) & 
                                  is.na(cnt_age_est_max))) %>% 
  mutate(age_range_mean = (cnt_age_est_min + cnt_age_est_max)/ 2, 
         age_group = ifelse(is.na(cnt_age_exact), 
                            ifelse(cnt_age_est_min < 75, 5 * floor(age_range_mean/5), 75), 
                            ifelse(cnt_age_exact < 75, 5 * floor(cnt_age_exact/5), 75))) %>% 
  select(part_id, age_group, everything(), -cnt_age_exact, 
         -cnt_age_est_min, -cnt_age_est_max, -age_range_mean,
         -cont_id, -cnt_gender, -frequency_multi, -phys_contact, -duration_multi) %>% 
  left_join(participants, by = "part_id") %>% # Adding in the age of the participant
  group_by(part_age, age_group) %>% # Grouping by participant age and age of person contacted
  summarise(total = sum(cnt_work)) %>% # Summing over all contacts
  left_join(number_people_df, by = c("part_age" = "age_group")) %>% # adding in the number of people surveyed in each age group
  mutate(per_cap = total/n) %>% # dividing total contacts by number surveyed in each age group
  select(part_age, age_group, per_cap) %>%
  spread(age_group, per_cap) # creating the matrix 
ages <- work_contacts$part_age 
work_contacts <- as.matrix(work_contacts[, -1]) # remove column of age group labels
work_contacts[is.na(work_contacts)] <- 0 # change the single NA at element [16, 14] which is an NA as 0 contacts between those 2 age groups were reported
row.names(work_contacts) <- ages
heatmap(work_contacts, Rowv = NA, Colv = NA)
saveRDS(work_contacts, file = "scripts/work_contact_matrix.rds")

# Generate Other Matrix
other_contacts <- polymod$contacts %>%
  filter(part_id %in% UK_ids, !(is.na(cnt_age_exact) & 
                                  !is.na(cnt_age_est_min) & 
                                  is.na(cnt_age_est_max))) %>% 
  mutate(age_range_mean = (cnt_age_est_min + cnt_age_est_max)/ 2, 
         age_group = ifelse(is.na(cnt_age_exact), 
                            ifelse(cnt_age_est_min < 75, 5 * floor(age_range_mean/5), 75), 
                            ifelse(cnt_age_exact < 75, 5 * floor(cnt_age_exact/5), 75))) %>% 
  select(part_id, age_group, everything(), -cnt_age_exact, 
         -cnt_age_est_min, -cnt_age_est_max, -age_range_mean,
         -cont_id, -cnt_gender, -frequency_multi, -phys_contact, -duration_multi) %>% 
  left_join(participants, by = "part_id") %>% # Adding in the age of the participant
  group_by(part_age, age_group) %>% # Grouping by participant age and age of person contacted
  summarise(total = sum(cnt_transport + cnt_leisure + cnt_otherplace)) %>% # Summing over all contacts
  left_join(number_people_df, by = c("part_age" = "age_group")) %>% # adding in the number of people surveyed in each age group
  mutate(per_cap = total/n) %>% # dividing total contacts by number surveyed in each age group
  select(part_age, age_group, per_cap) %>%
  spread(age_group, per_cap) # creating the matrix 
ages <- other_contacts$part_age 
other_contacts <- as.matrix(other_contacts[, -1]) # remove column of age group labels
other_contacts[is.na(other_contacts)] <- 0 # change the single NA at element [16, 14] which is an NA as 0 contacts between those 2 age groups were reported
row.names(other_contacts) <- ages
saveRDS(other_contacts, file = "scripts/other_contact_matrix.rds")

