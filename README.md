# Opentrons

## How to run this protocol

This repository includes two Opentrons Flex protocols:

- `ADH_assay/PQQ-ADH_assay.py`
- `ADH_assay/DCPIP_titrations.py`

### Run on the robot (Opentrons App)

1. Open the Opentrons App and connect to your Flex.
2. Import one of the protocol `.py` files above.
3. Confirm required labware/modules match your deck setup.
4. For `PQQ-ADH_assay.py`, set the runtime parameter `USE_TEMP_MODULE` as needed.
5. Start the run from the App.

### Optional: simulate before running

From the repository root:

```bash
opentrons_simulate ADH_assay/PQQ-ADH_assay.py
opentrons_simulate ADH_assay/DCPIP_titrations.py
```
