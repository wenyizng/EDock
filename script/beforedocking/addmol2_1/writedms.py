import chimera
from chimera import runCommand, openModels, MSMSModel
from WriteDMS import writeDMS
runCommand("surf")
surf = openModels.list(modelTypes=[MSMSModel])[0]
writeDMS(surf,"rec.ms")
