# deploy_model.py

import os
import torch
from mattersim.forcefield.potential import Potential
from mattersim.jit_compile_tools.jit_compile import deploy

CHECKPOINT_PATH = "./mattersim-v1.0.0-5M.pth"
DEPLOYED_MODEL_NAME = "deployed_mattersim.pth"
DEVICE = "cpu"

potential = Potential.from_checkpoint("./mattersim-v1.0.0-5M.pth", load_training_state=False)


deploy(
    potential,
    is_m3gnet_pretrained=True,
    deployed_model_name=DEPLOYED_MODEL_NAME,
    device=DEVICE,
)

print("Model has been successfully deployed!")