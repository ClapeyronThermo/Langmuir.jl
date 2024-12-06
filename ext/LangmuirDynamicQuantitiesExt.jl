module LangmuirDynamicQuantitiesExt

import Langmuir
import DynamicQuantities
using DynamicQuantities: @u_str

function Langmuir.Rgas(model::Langmuir.IsothermModel{<:DynamicQuantities.Quantity})
    return 8.31446261815324u"J/mol/K"
end

end #module