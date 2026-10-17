# Manual ICA review metadata

`frozen_ICA_restored_components.csv` records the component restorations encoded in `STEP7C_VERIFY_FINAL_ICA_DECISIONS.m` for the 14 participants that exceeded the >20% automatic ICA-rejection review trigger.

These are public dataset participant identifiers. The final rejected-component set for each reviewed participant is reconstructed as the automatic rejection set minus the listed restored components.

The interactive review script (`STEP7B_manual_review_flagged_ICA.m`) and the verification script (`STEP7C_VERIFY_FINAL_ICA_DECISIONS.m`) remain the authoritative implementation.
