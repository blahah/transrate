module Transrate

	class Writer

		require 'csv'

		def self.write name, data
			CSV.open(name, 'wb') do |csv|
				csv << ["metric", "value"]
				data.each_pair do |k, v|
					csv << [k, v]
				end
			end
		end

	end # Writer

end # Transrate
