library(tidyverse)

## Code from Hansen & Bowers (2009)
gotv <- haven::read_dta("NHrep_household.dta")
ivdl <- haven::read_dta("NHrep_individual.dta")

gotv <- gotv %>%
  rename_with(tolower) %>%
  filter(persons == 1, mailgrp == 0) %>%  
  mutate(
    v98_1 = recode(v98_1, `99` = 0),
    age1 = ifelse(age1miss == 1, NA, age1),
    age2 = ifelse(age2miss == 1, NA, age2),
    phongotv = ifelse(blood != 0, 1, phongotv)
  ) %>%
  rename(
    inperson_rand = persngrp,
    phone_rand = phongotv,
    turnout_98 = v98_1,
    age = age1,
    maj_party = majpty1,
    turnout_96 = v96_1_1,
    id = id1
    )

ivdl <- ivdl %>%
  filter(persons == 1) %>%  
  select(id1, cntany, pcntany) %>%
  rename(
    inperson = cntany, phone = pcntany, id = id1
  )

gotv <- left_join(gotv, ivdl, by = "id") %>%
  select(id, ward, turnout_98, inperson, phone, inperson_rand, phone_rand,
         age, maj_party, turnout_96)

gotv <- as.data.frame(gotv)
row.names(gotv) <- as.character(gotv$id)
newhaven <- gotv %>% select(-id)


save(newhaven, file = "../data/newhaven.rda", version = 2, compress = "bzip2")
